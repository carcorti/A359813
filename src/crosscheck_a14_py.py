#!/usr/bin/env python3
"""
A359813 — Run 3 cross-check indipendente
==========================================

Validazione incrociata dell'enumeratore del codice C `A359813.c` tramite
una implementazione Python completamente indipendente che usa `sympy.isprime`
(primalita' BPSW, deterministica fino a 2^64).

Scopo paper-grade: dimostrare che un prefisso completo dell'enumerazione
strutturale di `a(14)` produce lo stesso conteggio con:
  - `A359813` (C17, Miller-Rabin 12 basi, Montgomery form)
  - questo script (Python 3, sympy.isprime BPSW)

Prefisso scelto: pid = 75 (cifre 8, 0, 0).
Motivazione: sum_mod3 iniziale = 2, copre tutte le mask mod-3 nei 5^10 stati.

Uso:
    python3 crosscheck_a14_py.py                    # default pid=75
    python3 crosscheck_a14_py.py 75                 # esplicito
    python3 crosscheck_a14_py.py 75 --no-filters    # senza filtri modulari
    python3 crosscheck_a14_py.py --all              # tutti i 100 prefissi (lento!)

Output:
    - conteggio primi nel prefisso (deve corrispondere al C)
    - tempo di esecuzione
    - log su stdout

Firma: Sofia (Claude Opus 4.7), 19 aprile 2026.
AI-assisted source generated within the A359813 project workflow.
"""

import sys
import time
from itertools import product
from sympy import isprime


# ----- Costanti coerenti col codice C -----
SMALL_P = (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
LAST_DIGITS = (1, 3, 7, 9)
MIDDLE_DIGITS = (0, 2, 4, 6, 8)
FIRST_DIGITS = (2, 4, 6, 8)

# Mask mod 3: bit i corrisponde a LAST_DIGITS[i]
# mod3=0: last non puo' avere residuo 0 mod 3 (3,9) -> bit per 1,7 → 0b0101 = 5
# mod3=1: last non puo' avere residuo 2 mod 3 (nessuno ha residuo 2) → tutti → 0b1111 = 15
# mod3=2: last non puo' avere residuo 1 mod 3 (1,7) -> bit per 3,9 → 0b1010 = 10
LAST_OK_MASK = {0: 0b0101, 1: 0b1111, 2: 0b1010}


def decode_prefix3(pid):
    """Decodifica prefix_id in (pre[0], pre[1], pre[2]) — stessa formula del C."""
    assert 0 <= pid < 100, f"pid must be in [0,99], got {pid}"
    a = pid // 25
    r = pid % 25
    b = r // 5
    c = r % 5
    return (2 * (a + 1), 2 * b, 2 * c)


def rejected_by_small_filters(n):
    """Riproduce rejected_by_small_filters del C: include la guardia `n > p`."""
    for p in SMALL_P:
        if n > p and n % p == 0:
            return True
    return False


def count_prefix_a14(pid, use_filters=True, verbose=True):
    """
    Conta i primi del prefisso `pid` per a(14).
    
    Struttura: pre[0..2] + middle[0..9] + last, totale 14 cifre.
    middle_len = k - 4 = 10 per k=14.
    Stati middle = 5^10 = 9,765,625.
    """
    k = 14
    middle_len = k - 4  # = 10
    states = 5 ** middle_len
    
    pre = decode_prefix3(pid)
    prefix_sum_mod3 = sum(pre) % 3
    prefix_value = pre[0] * 10**(k-1) + pre[1] * 10**(k-2) + pre[2] * 10**(k-3)
    
    if verbose:
        print(f"  k = {k}")
        print(f"  prefix id = {pid}")
        print(f"  prefix digits = {pre}")
        print(f"  prefix value = {prefix_value:,}")
        print(f"  prefix sum_mod3 = {prefix_sum_mod3}")
        print(f"  middle_len = {middle_len}")
        print(f"  middle states = {states:,} (= 5^{middle_len})")
        print(f"  use_filters = {use_filters}")
        print()
    
    count = 0
    primes_checked = 0
    rejected_mod3 = 0
    rejected_small = 0
    
    # Enumerazione esplicita di tutti i 5^10 stati del middle, tramite itertools.
    t0 = time.time()
    progress_step = max(1, states // 20)
    
    for state_idx, middle in enumerate(product(MIDDLE_DIGITS, repeat=middle_len)):
        # Calcola sum_mod3 del candidato (prefix + middle + last ≡ prefix + middle mod 3, 
        # last aggiunta separatamente)
        middle_sum = sum(middle)
        sum_mod3 = (prefix_sum_mod3 + middle_sum) % 3
        mask = LAST_OK_MASK[sum_mod3]
        
        # Costruisce il "base_no_last" come prefix + middle (14 cifre meno quella finale)
        # middle occupa posizioni k-4..1 (posizione 0 riservata al last)
        # posizioni weight: middle[0] -> 10^(k-4)=10^10, middle[1] -> 10^9, ..., middle[9] -> 10^1
        base_no_last = prefix_value
        for t_idx, d in enumerate(middle):
            pos = k - 4 - t_idx  # k-4 down to 1
            base_no_last += d * 10**pos
        
        for li, last in enumerate(LAST_DIGITS):
            if ((mask >> li) & 1) == 0:
                rejected_mod3 += 1
                continue
            n = base_no_last + last
            if use_filters and rejected_by_small_filters(n):
                rejected_small += 1
                continue
            primes_checked += 1
            if isprime(n):
                count += 1
        
        if verbose and state_idx and state_idx % progress_step == 0:
            elapsed = time.time() - t0
            pct = 100 * state_idx / states
            eta = elapsed * (states - state_idx) / state_idx
            print(f"  [{pct:5.1f}%] state {state_idx:,}/{states:,} | "
                  f"found={count} | primes_tested={primes_checked:,} | "
                  f"elapsed={elapsed:.1f}s | eta={eta:.1f}s")
    
    elapsed = time.time() - t0
    
    return {
        'pid': pid,
        'pre': pre,
        'count': count,
        'primes_checked': primes_checked,
        'rejected_mod3': rejected_mod3,
        'rejected_small': rejected_small,
        'elapsed': elapsed,
    }


def main():
    # Parsing argomenti
    pid = 75
    use_filters = True
    run_all = False
    
    args = sys.argv[1:]
    i = 0
    while i < len(args):
        a = args[i]
        if a == '--no-filters':
            use_filters = False
        elif a == '--all':
            run_all = True
        elif a.isdigit():
            pid = int(a)
        else:
            print(f"Argomento sconosciuto: {a}", file=sys.stderr)
            print(__doc__, file=sys.stderr)
            return 2
        i += 1
    
    print("=" * 68)
    print("A359813 — Run 3 cross-check indipendente (Python + sympy.isprime)")
    print("=" * 68)
    print(f"Target: un prefisso di a(14) completo")
    print()
    
    if run_all:
        print("Modalita' --all: lancio tutti i 100 prefissi (molto lento!)")
        total = 0
        t0 = time.time()
        for pid in range(100):
            r = count_prefix_a14(pid, use_filters=use_filters, verbose=False)
            total += r['count']
            print(f"  prefix {pid:3d} (pre={r['pre']}): {r['count']:,} primi in {r['elapsed']:.1f}s")
        elapsed = time.time() - t0
        print()
        print(f"TOTALE a(14) layer = {total:,}")
        print(f"ORACLE a(14) layer = 310498347")
        print(f"Match: {'YES' if total == 310498347 else 'NO'}")
        print(f"Elapsed totale: {elapsed:.1f}s = {elapsed/60:.1f}min")
    else:
        r = count_prefix_a14(pid, use_filters=use_filters, verbose=True)
        print()
        print("=" * 68)
        print("RISULTATO FINALE")
        print("=" * 68)
        print(f"  prefix id        : {r['pid']}")
        print(f"  prefix digits    : {r['pre']}")
        print(f"  primi trovati    : {r['count']:,}")
        print(f"  primi testati    : {r['primes_checked']:,}")
        print(f"  rigettati mod 3  : {r['rejected_mod3']:,}")
        print(f"  rigettati small  : {r['rejected_small']:,}")
        print(f"  elapsed          : {r['elapsed']:.1f}s = {r['elapsed']/60:.1f}min")
        print()
        print("Per confrontare col C, modificare count_layer_k nel sorgente per")
        print(f"iterare solo su pid == {r['pid']} (oppure aggiungere flag --only-prefix).")
        print()
        print("In alternativa, lanciare --all per ottenere a(14) completo (310,498,347).")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
