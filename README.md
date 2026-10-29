# Computational archive

for *On the Lisoněk–Salt conjecture: tilings and Kochen–Specker pairs from Hadamard matrices* (Vijaykumar Singh).

`code/COMPUTATIONS.md` lists every computational statement of the paper with the script that produces it, its
inputs and normalizations, the number of cases and the outcome. `REPRODUCTION.md` records a complete rerun from this
archive. The computations confirm statements proved in the paper, settle the finite cases that the paper marks as
computational (Theorem 12.1, and the case lambda = 2, p = 7 in Remark 11.7), and support Remark 8.14.

## Contents

| Path | |
|---|---|
| `code/*.py`, `code/quadric/*.py` | Python scripts; run from `code/` (or `code/quadric/` where the manifest says so) |
| `code/h24_min.c` | C program for the order-24 search (Theorem 9.6): `gcc -O2 -o h24_min h24_min.c`, then `python3 h24_min.py h24_rows.txt && ./h24_min 4 < h24_rows.txt` |
| `code/quadric/dimcheck_cluster.py` | sharded version of the q = 25 dimension-4 check (Corollary 8.13); `selftest` runs q = 9 |
| `code/lean/Conic.lean` | Lean 4 proof of the nonvanishing used in Theorem 6.1 (Lean v4.34.1, Mathlib) |
| `code/sloane/`, `code/butson/` | input matrices: Sloane's Hadamard matrices of order 16, Butson matrices BH(10,5), BH(14,7) |
| `code/witnesses/` | explicit partitions for Table 2 |
| `code/logs/`, `code/quadric/logs/`, `code/*.log` | logs of the runs cited in the manifest, including the long ones |

## Requirements

Python 3 (tested with 3.14), `python-sat` (PySAT; tested with 1.9.dev17, CaDiCaL 1.5.3) and `numpy` (tested with
2.5); a C compiler for `h24_min.c`; Lean 4 with Mathlib for `Conic.lean`. Random choices use fixed seeds stated in
the scripts and the manifest.

## Licence

Except as noted below, the contents of this archive are dedicated to the public domain under CC0 1.0 Universal (see
`LICENSE`). The matrices in `code/sloane/` are from N. J. A. Sloane's library of Hadamard matrices
(http://neilsloane.com/hadamard/), and those in `code/butson/` are from the tables of Butson Hadamard matrices of
P. Lampio, P. R. J. Östergård and F. Szöllősi (https://wiki.aalto.fi/display/Butson/Matrices+up+to+monomial+equivalence); they are included for reproducibility and keep their original
terms.
