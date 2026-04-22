# OEIS A359813 — Computation of New Terms

This repository contains the source code, data files, selected logs, and paper associated with the computation of two new terms of [OEIS A359813](https://oeis.org/A359813).

> **A359813** — Number of primes < 10^n with exactly one odd decimal digit.

## Main results

|  n |                a(n) | Status         |
| -: | ------------------: | -------------- |
| 17 |  **40 386 401 580** | New, this work |
| 18 | **190 186 145 230** | New, this work |

Terms `a(1)` through `a(16)` (Zhining Yang, 2023) were independently re-verified.
The full sequence through `a(18)` is provided in [`data/b359813.txt`](data/b359813.txt).

## Repository contents

```text
A359813/
├── README.md
├── LICENSE
├── CITATION.cff
├── src/       # C source code, build files, and run scripts
├── data/      # sequence data
├── logs/      # selected execution and validation logs
└── paper/     # PDF and LaTeX source of the accompanying paper
```

The accompanying paper is available as [`paper/A359813.pdf`](paper/A359813.pdf) and [`paper/A359813.tex`](paper/A359813.tex).

## Mathematical structure of the search space

A prime below `10^n` with exactly one odd decimal digit must have that odd digit in the units position.
Indeed:

- if the last digit is even, the number is divisible by `2`;
- if the last digit is `5`, the number is divisible by `5`.

Hence the candidate set is

```text
{2,4,6,8} × {0,2,4,6,8}^(n−2) × {1,3,7,9}
```

with cardinality

```text
16 · 5^(n−2)
```

This reduces the search space exponentially relative to the full interval below `10^n`.

## Algorithmic overview

Each candidate is processed through the following pipeline:

1. **Small-prime filtering**  
   Modular sieving with primes `p ≤ 47`, using incrementally updated residues.

2. **Primality testing**  
   Deterministic Miller–Rabin with 12 bases
   `{2,3,5,7,11,13,17,19,23,29,31,37}`
   for the explored numerical range.

3. **Arithmetic kernels**
   - Montgomery modular multiplication (primary backend)
   - Independent `__uint128_t` backend (cross-check)

4. **Parallelisation**  
   OpenMP over the 100 three-digit prefixes.

5. **Checkpointing**  
   Binary checkpoints with CRC32 integrity checks, enabling safe resume.

The implementation is written in C17 and uses only libc plus OpenMP.
Full algorithmic details are documented in the accompanying paper.

## Hardware and software environment

- CPU: AMD Ryzen 9 7940HS (8 cores / 16 threads)
- RAM: 64 GB DDR5
- OS: Linux Mint 22.3 (kernel 6.14)
- Compiler: GCC 13.3.0 with `-O3 -march=native -fopenmp -std=c17`

## Performance

| Computation                          | Threads | Wall-time |
| ------------------------------------ | ------: | --------: |
| Oracle test (`a(1)`…`a(16)`)         |      16 |   ~43 min |
| `a(17)`, Montgomery backend          |      16 |  2 h 20 m |
| `a(17)`, `__uint128_t` backend       |      16 |  4 h 12 m |
| `a(18)`, Montgomery backend          |      14 | 12 h 49 m |

Wall-times are hardware-dependent.

## Validation

The reported results were validated through three independent layers.

### 1. Oracle verification

Recomputation of `a(1)` through `a(16)` matches the OEIS b-file exactly.

### 2. Backend independence

`a(17)` was computed with both:

- Montgomery arithmetic
- native `__uint128_t`

The two runs produced identical results.

### 3. Independent enumerator cross-check

Selected prefixes of `a(14)` were cross-checked with a separate Python enumerator using `sympy.isprime` / BPSW.

## Heuristic consistency checks

These are consistency indicators, not proofs.

- **Prime number theorem scaling**  
  `a(18)/a(17) ≈ 4.7092`, expected `≈ 4.7150`, deviation `≈ 0.12%`

- **Density heuristic**  
  Agreement with `π(10^n)/2^(n−1)` at about the 1% level.

## Reproducibility

The scripts in `src/` are intended for Linux environments with GCC and OpenMP support.

### Build

```bash
cd src/
make release
```

### Mandatory validation run

```bash
./oracle_test.sh
```

### Compute `a(17)`

```bash
./campaign_a17.sh
```

### Compute `a(18)`

```bash
./campaign_a18.sh
```

## Command-line interface

```text
Usage: A359813 [--oracle] [--run N] [--threads T] [--resume]
               [--backend mont|u128] [--paranoid12]
               [--checkpoint PATH] [--log PATH]
               [--quiet|--trace]
```

## License and citation

This project is released under the MIT License; see [`LICENSE`](LICENSE).
Citation metadata are provided in [`CITATION.cff`](CITATION.cff).

## Acknowledgements

This project was developed within an AI-assisted workflow.
All published numerical results were independently validated through:

- OEIS oracle verification
- backend independence
- independent cross-check computations

No published result relies on a single computational path.
