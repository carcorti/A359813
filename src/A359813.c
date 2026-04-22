/*
 * A359813.c — C17/OpenMP campaign code for OEIS A359813
 *
 * Number of primes < 10^n with exactly one odd decimal digit.
 *
 * Author:   Carlo Stradiotti (AI-assisted)
 * Platform: AMD Ryzen 9 7940HS, Linux Mint 22.3, GCC 13.3 + OpenMP
 * License:  MIT (see LICENSE in repository root)
 *
 * Primary result : a(17) = 40386401580    (validated bit-per-bit Run1=Run2)
 * Secondary      : a(18) = 190186145230   (single-run, Wilson rel = 0.75%)
 *
 * Verification: MR deterministic with 12 bases {2,3,5,7,11,13,17,19,23,29,31,37}
 * (deterministic up to n < 3.3e23, Jaeschke 1993 / OEIS A014233).
 *
 * Build:   make release     (or: make debug / profile / sanitize)
 * Oracle:  ./oracle_test.sh
 * Run 17:  ./campaign_a17.sh
 * Run 18:  ./A359813 --run 18 --threads 14 --backend mont --paranoid12 \
 *            --checkpoint A359813_state_a18.bin --log a18.log
 */

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define A359813_VERSION 1u
#define MAX_N 18u
#define PREFIX_DEPTH 3u
#define N_PREFIXES 100u
#define N_SMALLP 12u
#define MR7_NBASES 7u
#define MR12_NBASES 12u
#define LOG_QUIET 0
#define LOG_INFO  1
#define LOG_TRACE 2

static const uint32_t SMALL_P[N_SMALLP] = {7,11,13,17,19,23,29,31,37,41,43,47};
static const uint64_t MR_BASES_7[MR7_NBASES]  = {2,3,5,7,11,13,17};
static const uint64_t MR_BASES_12[MR12_NBASES] = {2,3,5,7,11,13,17,19,23,29,31,37};
static const uint64_t ORACLE_A[18] = {
    0ULL,
    3ULL, 12ULL, 45ULL, 171ULL, 619ULL, 2560ULL, 10774ULL, 46708ULL,
    202635ULL, 904603ULL, 4073767ULL, 18604618ULL, 85445767ULL,
    395944114ULL, 1837763447ULL, 8600149593ULL,
    40386401580ULL /* a(17), validated by Run1=Run2 bit-per-bit, Apr 19 2026 */
};

typedef struct {
    uint32_t target_k;
    uint32_t n_prefixes;
    uint8_t  completed[(N_PREFIXES + 7u) / 8u];
    uint64_t partial_counts[N_PREFIXES];
    uint32_t crc32;
} checkpoint_t;

typedef struct {
    uint64_t pow10[MAX_N + 1u];
    uint64_t pow5[MAX_N + 1u];
    uint32_t contrib[N_SMALLP][MAX_N][5];
    uint32_t last_add[N_SMALLP][4];
    uint8_t  last_ok_mask[3];
} tables_t;

typedef struct {
    uint64_t n;
    uint64_t nprime;
    uint64_t r2;
} mont_ctx_t;

typedef struct {
    unsigned target_n;
    unsigned threads;
    bool resume;
    bool use_mont;
    bool paranoid_12bases;
    int log_level;
    const char *checkpoint_path;
    const char *log_path;
} run_config_t;

static FILE *g_log = NULL;

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void log_msg(int level, const run_config_t *cfg, const char *fmt, ...) {
    if (cfg->log_level < level) return;
    va_list ap;
    va_start(ap, fmt);
    vfprintf(stdout, fmt, ap);
    fflush(stdout);
    va_end(ap);
    if (g_log) {
        va_start(ap, fmt);
        vfprintf(g_log, fmt, ap);
        fflush(g_log);
        va_end(ap);
    }
}

static void tables_init(tables_t *tb) {
    memset(tb, 0, sizeof(*tb));
    tb->pow10[0] = 1ULL;
    tb->pow5[0] = 1ULL;
    for (unsigned i = 1; i <= MAX_N; ++i) {
        tb->pow10[i] = tb->pow10[i - 1] * 10ULL;
        tb->pow5[i] = tb->pow5[i - 1] * 5ULL;
    }
    tb->last_ok_mask[0] = 0x5u; /* {1,7} */
    tb->last_ok_mask[1] = 0xFu; /* {1,3,7,9} */
    tb->last_ok_mask[2] = 0xAu; /* {3,9} */

    for (unsigned pi = 0; pi < N_SMALLP; ++pi) {
        uint32_t p = SMALL_P[pi];
        for (unsigned pos = 0; pos <= MAX_N - 1u; ++pos) {
            uint32_t w = (uint32_t)(tb->pow10[pos] % p);
            for (unsigned idx = 0; idx < 5; ++idx) {
                tb->contrib[pi][pos][idx] = (uint32_t)((2u * idx * w) % p);
            }
        }
        tb->last_add[pi][0] = 1u % p;
        tb->last_add[pi][1] = 3u % p;
        tb->last_add[pi][2] = 7u % p;
        tb->last_add[pi][3] = 9u % p;
    }
}

static void decode_prefix3(uint32_t prefix_id, uint8_t out[3]) {
    uint32_t a = prefix_id / 25u;
    uint32_t r = prefix_id % 25u;
    uint32_t b = r / 5u;
    uint32_t c = r % 5u;
    out[0] = (uint8_t)(2u * (a + 1u));
    out[1] = (uint8_t)(2u * b);
    out[2] = (uint8_t)(2u * c);
}

static inline uint32_t add_mod_u32(uint32_t a, uint32_t b, uint32_t m) {
    uint32_t s = a + b;
    return (s >= m) ? (s - m) : s;
}

static inline uint32_t sub_mod_u32(uint32_t a, uint32_t b, uint32_t m) {
    return (a >= b) ? (a - b) : (a + m - b);
}

static inline uint64_t mulmod_u128(uint64_t a, uint64_t b, uint64_t mod) {
    return (uint64_t)(((__uint128_t)a * (__uint128_t)b) % mod);
}

static inline uint64_t powmod_u128(uint64_t a, uint64_t e, uint64_t mod) {
    uint64_t r = 1 % mod;
    while (e) {
        if (e & 1ULL) r = mulmod_u128(r, a, mod);
        a = mulmod_u128(a, a, mod);
        e >>= 1ULL;
    }
    return r;
}

static bool mont_ctx_init(mont_ctx_t *ctx, uint64_t n) {
    if ((n & 1ULL) == 0ULL) return false;
    ctx->n = n;
    uint64_t x = 1ULL;
    for (unsigned i = 0; i < 6; ++i) x *= 2ULL - n * x;
    ctx->nprime = (uint64_t)(0ULL - x);
    uint64_t r_mod_n = ((uint64_t)0ULL - n) % n; /* 2^64 mod n */
    ctx->r2 = mulmod_u128(r_mod_n, r_mod_n, n);
    return true;
}

static inline uint64_t mont_reduce(__uint128_t t, const mont_ctx_t *ctx) {
    uint64_t m = (uint64_t)t * ctx->nprime;
    __uint128_t u = (t + (__uint128_t)m * ctx->n) >> 64;
    uint64_t r = (uint64_t)u;
    return (r >= ctx->n) ? (r - ctx->n) : r;
}

static inline uint64_t mont_mul(uint64_t a, uint64_t b, const mont_ctx_t *ctx) {
    return mont_reduce((__uint128_t)a * b, ctx);
}

static inline uint64_t mont_encode(uint64_t x, const mont_ctx_t *ctx) {
    return mont_reduce((__uint128_t)x * ctx->r2, ctx);
}

static inline uint64_t mont_decode(uint64_t x, const mont_ctx_t *ctx) {
    return mont_reduce((__uint128_t)x, ctx);
}

static uint64_t mont_pow(uint64_t a, uint64_t e, const mont_ctx_t *ctx) {
    uint64_t one = mont_encode(1ULL, ctx);
    uint64_t base = mont_encode(a % ctx->n, ctx);
    uint64_t r = one;
    while (e) {
        if (e & 1ULL) r = mont_mul(r, base, ctx);
        base = mont_mul(base, base, ctx);
        e >>= 1ULL;
    }
    return r;
}

static bool mr_core_u128(uint64_t n, const uint64_t *bases, unsigned nbases) {
    uint64_t d = n - 1ULL;
    unsigned s = 0;
    while ((d & 1ULL) == 0ULL) { d >>= 1ULL; ++s; }
    for (unsigned i = 0; i < nbases; ++i) {
        uint64_t a = bases[i] % n;
        if (a == 0ULL) continue;
        uint64_t x = powmod_u128(a, d, n);
        if (x == 1ULL || x == n - 1ULL) continue;
        bool ok = false;
        for (unsigned r = 1; r < s; ++r) {
            x = mulmod_u128(x, x, n);
            if (x == n - 1ULL) { ok = true; break; }
        }
        if (!ok) return false;
    }
    return true;
}

static bool mr_core_mont(uint64_t n, const uint64_t *bases, unsigned nbases) {
    mont_ctx_t ctx;
    if (!mont_ctx_init(&ctx, n)) return false;
    uint64_t d = n - 1ULL;
    unsigned s = 0;
    while ((d & 1ULL) == 0ULL) { d >>= 1ULL; ++s; }
    uint64_t oneM = mont_encode(1ULL, &ctx);
    uint64_t nm1M = mont_encode(n - 1ULL, &ctx);
    for (unsigned i = 0; i < nbases; ++i) {
        uint64_t a = bases[i] % n;
        if (a == 0ULL) continue;
        uint64_t x = mont_pow(a, d, &ctx);
        if (x == oneM || x == nm1M) continue;
        bool ok = false;
        for (unsigned r = 1; r < s; ++r) {
            x = mont_mul(x, x, &ctx);
            if (x == nm1M) { ok = true; break; }
        }
        if (!ok) return false;
    }
    return true;
}

static bool is_prime_u64(uint64_t n, bool use_mont, bool paranoid12) {
    static const uint32_t SMALL0[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    if (n < 2ULL) return false;
    for (size_t i = 0; i < sizeof(SMALL0)/sizeof(SMALL0[0]); ++i) {
        if (n == SMALL0[i]) return true;
        if (n % SMALL0[i] == 0ULL) return false;
    }
    const uint64_t *bases = paranoid12 ? MR_BASES_12 : MR_BASES_7;
    unsigned nbases = paranoid12 ? MR12_NBASES : MR7_NBASES;
    return use_mont ? mr_core_mont(n, bases, nbases) : mr_core_u128(n, bases, nbases);
}

static uint32_t crc32_ieee(const void *data, size_t len) {
    static uint32_t table[256];
    static bool inited = false;
    if (!inited) {
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (unsigned j = 0; j < 8; ++j) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            table[i] = c;
        }
        inited = true;
    }
    uint32_t c = 0xFFFFFFFFu;
    const uint8_t *p = (const uint8_t *)data;
    for (size_t i = 0; i < len; ++i) c = table[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

static void checkpoint_zero(checkpoint_t *ckp, uint32_t target_k) {
    memset(ckp, 0, sizeof(*ckp));
    ckp->target_k = target_k;
    ckp->n_prefixes = N_PREFIXES;
}

static inline bool checkpoint_done(const checkpoint_t *ckp, uint32_t pid) {
    return (ckp->completed[pid >> 3u] >> (pid & 7u)) & 1u;
}

static inline void checkpoint_mark_done(checkpoint_t *ckp, uint32_t pid) {
    ckp->completed[pid >> 3u] |= (uint8_t)(1u << (pid & 7u));
}

static bool checkpoint_save_atomic(const char *path, checkpoint_t *ckp) {
    char tmp[1024];
    snprintf(tmp, sizeof(tmp), "%s.tmp", path);
    ckp->crc32 = 0u;
    ckp->crc32 = crc32_ieee(ckp, sizeof(*ckp) - sizeof(uint32_t));
    FILE *f = fopen(tmp, "wb");
    if (!f) return false;
    if (fwrite(ckp, sizeof(*ckp), 1, f) != 1) { fclose(f); return false; }
    fflush(f);
    fclose(f);
    return rename(tmp, path) == 0;
}

static bool checkpoint_load(const char *path, checkpoint_t *ckp) {
    FILE *f = fopen(path, "rb");
    if (!f) return false;
    if (fread(ckp, sizeof(*ckp), 1, f) != 1) { fclose(f); return false; }
    fclose(f);
    uint32_t got = ckp->crc32;
    ckp->crc32 = 0u;
    uint32_t want = crc32_ieee(ckp, sizeof(*ckp) - sizeof(uint32_t));
    ckp->crc32 = got;
    return got == want && ckp->n_prefixes == N_PREFIXES;
}

static inline uint32_t last_digit_by_idx(unsigned idx) {
    static const uint32_t d[4] = {1u,3u,7u,9u};
    return d[idx];
}

static inline bool rejected_by_small_filters(uint64_t c, const uint32_t residues[N_SMALLP], unsigned last_idx, const tables_t *tb) {
    for (unsigned pi = 0; pi < N_SMALLP; ++pi) {
        uint32_t p = SMALL_P[pi];
        uint32_t r = add_mod_u32(residues[pi], tb->last_add[pi][last_idx], p);
        if (c > (uint64_t)p && r == 0u) return true;
    }
    return false;
}

static uint64_t count_layer_small(unsigned k, const tables_t *tb, const run_config_t *cfg) {
    (void)tb;
    uint64_t count = 0;
    if (k == 1u) return 3ULL;
    if (k == 2u) {
        const uint32_t firsts[4] = {2,4,6,8};
        const uint32_t lasts[4] = {1,3,7,9};
        for (unsigned i = 0; i < 4; ++i) for (unsigned j = 0; j < 4; ++j) {
            uint64_t c = 10ULL * firsts[i] + lasts[j];
            if (is_prime_u64(c, cfg->use_mont, cfg->paranoid_12bases)) ++count;
        }
        return count;
    }
    if (k == 3u) {
        const uint32_t firsts[4] = {2,4,6,8};
        const uint32_t mids[5] = {0,2,4,6,8};
        const uint32_t lasts[4] = {1,3,7,9};
        for (unsigned i = 0; i < 4; ++i) for (unsigned m = 0; m < 5; ++m) for (unsigned j = 0; j < 4; ++j) {
            uint64_t c = 100ULL * firsts[i] + 10ULL * mids[m] + lasts[j];
            if (is_prime_u64(c, cfg->use_mont, cfg->paranoid_12bases)) ++count;
        }
        return count;
    }
    return 0ULL;
}

static uint64_t count_layer_k(unsigned k, const tables_t *tb, const run_config_t *cfg, checkpoint_t *ckp) {
    if (k <= 3u) return count_layer_small(k, tb, cfg);
    const unsigned middle_len = k - 4u;
    volatile int prefixes_done = 0;
    double t0 = now_sec();
    uint64_t total = 0ULL;

#pragma omp parallel for schedule(dynamic,1) reduction(+:total)
    for (uint32_t pid = 0; pid < N_PREFIXES; ++pid) {
        if (ckp && checkpoint_done(ckp, pid)) {
            total += ckp->partial_counts[pid];
            continue;
        }

        uint8_t pre[3];
        decode_prefix3(pid, pre);
        uint64_t prefix_value = (uint64_t)pre[0] * tb->pow10[k - 1u]
                              + (uint64_t)pre[1] * tb->pow10[k - 2u]
                              + (uint64_t)pre[2] * tb->pow10[k - 3u];
        uint32_t prefix_res[N_SMALLP];
        for (unsigned pi = 0; pi < N_SMALLP; ++pi) {
            uint32_t p = SMALL_P[pi];
            prefix_res[pi] = 0u;
            prefix_res[pi] = add_mod_u32(prefix_res[pi], (uint32_t)((pre[0] * (tb->pow10[k - 1u] % p)) % p), p);
            prefix_res[pi] = add_mod_u32(prefix_res[pi], (uint32_t)((pre[1] * (tb->pow10[k - 2u] % p)) % p), p);
            prefix_res[pi] = add_mod_u32(prefix_res[pi], (uint32_t)((pre[2] * (tb->pow10[k - 3u] % p)) % p), p);
        }
        uint8_t idx[16] = {0};
        uint32_t cur_res[N_SMALLP];
        memcpy(cur_res, prefix_res, sizeof(cur_res));
        uint64_t middle_value = 0ULL;
        unsigned sum_mod3 = (pre[0] + pre[1] + pre[2]) % 3u;
        uint64_t local = 0ULL;
        const uint64_t states = tb->pow5[middle_len];

        for (uint64_t state = 0; state < states; ++state) {
            uint64_t base_no_last = prefix_value + middle_value;
            uint8_t mask = tb->last_ok_mask[sum_mod3];
            assert(mask != 0u);
            for (unsigned li = 0; li < 4; ++li) {
                if (((mask >> li) & 1u) == 0u) continue;
                uint64_t c = base_no_last + last_digit_by_idx(li);
                if (rejected_by_small_filters(c, cur_res, li, tb)) continue;
                if (is_prime_u64(c, cfg->use_mont, cfg->paranoid_12bases)) ++local;
            }
            if (state + 1u == states) break;
            for (int t = (int)middle_len - 1; t >= 0; --t) {
                unsigned pos = 1u + (unsigned)(middle_len - 1 - (unsigned)t);
                uint64_t weight = tb->pow10[pos];
                if (idx[t] < 4u) {
                    idx[t]++;
                    middle_value += 2ULL * weight;
                    sum_mod3 += 2u;
                    if (sum_mod3 >= 3u) sum_mod3 -= 3u;
                    for (unsigned pi = 0; pi < N_SMALLP; ++pi) cur_res[pi] = add_mod_u32(cur_res[pi], tb->contrib[pi][pos][1], SMALL_P[pi]);
                    break;
                } else {
                    unsigned old = idx[t];
                    idx[t] = 0u;
                    middle_value -= (uint64_t)(2u * old) * weight;
                    unsigned delta = (2u * old) % 3u;
                    sum_mod3 = (sum_mod3 + 3u - delta) % 3u;
                    for (unsigned pi = 0; pi < N_SMALLP; ++pi) cur_res[pi] = sub_mod_u32(cur_res[pi], tb->contrib[pi][pos][old], SMALL_P[pi]);
                }
            }
        }

        total += local;
        if (ckp) {
#pragma omp critical(checkpoint_io)
            {
                ckp->partial_counts[pid] = local;
                checkpoint_mark_done(ckp, pid);
                checkpoint_save_atomic(cfg->checkpoint_path, ckp);
                prefixes_done++;
                if (cfg->log_level >= LOG_INFO) {
                    double elapsed = now_sec() - t0;
                    double pct = 100.0 * (double)prefixes_done / (double)N_PREFIXES;
                    double eta = (prefixes_done > 0) ? elapsed * (double)(N_PREFIXES - prefixes_done) / (double)prefixes_done : 0.0;
                    log_msg(LOG_INFO, cfg, "[k=%u] prefix %3u/%u done | %.1f%% | partial=%" PRIu64 " | elapsed=%.1fs | eta=%.1fs\n",
                            k, prefixes_done, N_PREFIXES, pct, local, elapsed, eta);
                }
            }
        }
    }
    return total;
}

static uint64_t count_upto_oracle(unsigned n, const tables_t *tb, const run_config_t *cfg) {
    uint64_t cum = 0ULL;
    for (unsigned k = 1; k <= n; ++k) {
        uint64_t layer = count_layer_k(k, tb, cfg, NULL);
        cum += layer;
        log_msg(LOG_INFO, cfg, "oracle: k=%u layer=%" PRIu64 " cumulative=%" PRIu64 "\n", k, layer, cum);
    }
    return cum;
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [--oracle] [--run N] [--threads T] [--resume] [--backend mont|u128] [--paranoid12] [--checkpoint PATH] [--log PATH] [--quiet|--trace]\n",
        prog);
}

int main(int argc, char **argv) {
    run_config_t cfg = {
        .target_n = 17u,
        .threads = 16u,
        .resume = false,
        .use_mont = true,
        .paranoid_12bases = false,
        .log_level = LOG_INFO,
        .checkpoint_path = "A359813_state.bin",
        .log_path = "a17.log"
    };
    bool do_oracle = false;
    bool do_run = false;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--oracle") == 0) do_oracle = true;
        else if (strcmp(argv[i], "--run") == 0 && i + 1 < argc) { cfg.target_n = (unsigned)strtoul(argv[++i], NULL, 10); do_run = true; }
        else if (strcmp(argv[i], "--threads") == 0 && i + 1 < argc) cfg.threads = (unsigned)strtoul(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--resume") == 0) cfg.resume = true;
        else if (strcmp(argv[i], "--backend") == 0 && i + 1 < argc) {
            ++i;
            if (strcmp(argv[i], "mont") == 0) cfg.use_mont = true;
            else if (strcmp(argv[i], "u128") == 0) cfg.use_mont = false;
            else { usage(argv[0]); return 2; }
        } else if (strcmp(argv[i], "--paranoid12") == 0) cfg.paranoid_12bases = true;
        else if (strcmp(argv[i], "--checkpoint") == 0 && i + 1 < argc) cfg.checkpoint_path = argv[++i];
        else if (strcmp(argv[i], "--log") == 0 && i + 1 < argc) cfg.log_path = argv[++i];
        else if (strcmp(argv[i], "--quiet") == 0) cfg.log_level = LOG_QUIET;
        else if (strcmp(argv[i], "--trace") == 0) cfg.log_level = LOG_TRACE;
        else { usage(argv[0]); return 2; }
    }
    if (!do_oracle && !do_run) { usage(argv[0]); return 2; }
    if (cfg.target_n < 1u || cfg.target_n > MAX_N) {
        fprintf(stderr, "target_n must be in [1,%u]\n", MAX_N);
        return 2;
    }

    if (cfg.log_path) g_log = fopen(cfg.log_path, do_oracle ? "w" : "a");
    omp_set_num_threads((int)cfg.threads);

    tables_t tb;
    tables_init(&tb);

    if (do_oracle) {
        double t0 = now_sec();
        uint64_t got = 0ULL;
        for (unsigned n = 1; n <= 16u; ++n) {
            got = count_upto_oracle(n, &tb, &cfg);
            if (got != ORACLE_A[n]) {
                fprintf(stderr, "ORACLE FAIL at n=%u: got=%" PRIu64 " expected=%" PRIu64 "\n", n, got, ORACLE_A[n]);
                if (g_log) fclose(g_log);
                return 1;
            }
            log_msg(LOG_INFO, &cfg, "ORACLE OK n=%u value=%" PRIu64 "\n", n, got);
        }
        log_msg(LOG_INFO, &cfg, "oracle suite passed in %.3f s\n", now_sec() - t0);
        if (g_log) fclose(g_log);
        return 0;
    }

    if (do_run) {
        double t0 = now_sec();
        uint64_t result = 0ULL;
        if (cfg.target_n <= 16u) {
            result = ORACLE_A[cfg.target_n];
            log_msg(LOG_INFO, &cfg, "Using oracle value for a(%u) = %" PRIu64 "\n", cfg.target_n, result);
        } else {
            checkpoint_t ckp;
            bool loaded = false;
            if (cfg.resume) loaded = checkpoint_load(cfg.checkpoint_path, &ckp);
            if (!loaded) checkpoint_zero(&ckp, cfg.target_n);
            else if (ckp.target_k != cfg.target_n) checkpoint_zero(&ckp, cfg.target_n);
            uint64_t layer = count_layer_k(cfg.target_n, &tb, &cfg, &ckp);
            result = ORACLE_A[cfg.target_n - 1u] + layer; /* production logic: ORACLE_A[n-1] + layer_n */
            log_msg(LOG_INFO, &cfg, "layer delta(%u) = %" PRIu64 "\n", cfg.target_n, layer);
            if (cfg.target_n == 17u || cfg.target_n == 18u) {
                const long double wilson = (cfg.target_n == 17u)
                    ? (2623557157654233.0L / 65536.0L)        /* pi(10^17) / 2^16 */
                    : (24739954287740860.0L / 131072.0L);     /* pi(10^18) / 2^17 */
                long double rel = (result > wilson) ? ((long double)result - wilson) / (long double)result
                                                   : (wilson - (long double)result) / (long double)result;
                log_msg(LOG_INFO, &cfg, "Wilson sanity: a(%u)=%" PRIu64 " vs %.0Lf | rel=%.6Lf\n", cfg.target_n, result, wilson, rel);
                if (rel > 0.02L) {
                    fprintf(stderr, "FATAL: Wilson sanity failed (rel=%.6Lf > 0.02)\n", rel);
                    if (g_log) fclose(g_log);
                    return 3;
                }
            }
        }
        log_msg(LOG_INFO, &cfg, "RESULT a(%u) = %" PRIu64 " | elapsed=%.3f s | backend=%s | bases=%u\n",
                cfg.target_n, result, now_sec() - t0, cfg.use_mont ? "mont" : "u128", cfg.paranoid_12bases ? 12u : 7u);
        if (g_log) fclose(g_log);
        return 0;
    }

    if (g_log) fclose(g_log);
    return 0;
}
