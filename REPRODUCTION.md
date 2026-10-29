# Reproduction record

Rerun on 29 October 2026 from a clean copy of `code/` (Python 3.14.7, python-sat 1.9.dev17, numpy 2.5.3, Lean
v4.34.1 with Mathlib), each script with the arguments of `code/COMPUTATIONS.md`. Every outcome agrees with the
manifest.

| Manifest row | Command (from `code/` unless stated) | Outcome of the rerun |
|---|---|---|
| Thm 4.1, colourings of Thm 7.4 | `p3mod4_proof_check.py` | identities hold for the 8 primes p ≡ 3 (mod 4), p ≤ 47; colourings valid on every basis for p = 3, 7, 11, on 15 sampled hyperplanes for p = 19, 23, and through all values of β for every p ≡ 3 (mod 4) up to 79 (stopped at the 3-hour limit during p = 83, beyond the manifest's range) |
| Cor 5.2, Thm 6.1 | `verify.py` | p(2p−1) KS hyperplanes for p ≤ 13 (15, 45, 91, 231, 325); no subspace of dim 1 or 2 without a weight-1 slice |
| Prop 8.1, Thm 8.2 | `quadric/quadric_gh.py` | GH and rank 2m+2 for q = 3, 5, 9, 25, 27, 49 |
| Thm 8.5 | `quadric/quadric_weights.py` | weight classes and counts as stated |
| Thm 8.10 (m = 2) | `quadric/separable.py 3 2`, `5 2` | every non-KS hyperplane coloured; 0 needing another colouring |
| Remark after Prop 8.9 | `quadric/minbases_check.py 3 2`, `5 2` | 45,255 and 3,070,893 subspaces; 16 and 36 without a weight-1 slice, all colourable |
| Thm 8.11 | `quadric/e1_check.py p m 0` for q = 3, 5, 7, 9, 25 | 0 mismatches; KS of dim m+2: 6, 20, 42, 432, 12,000; 2,664 and 42,176 subspaces for q = 9, 25 |
| Thm 8.11, q = 27 | `quadric/e1_check.py 3 3 0` | 0 mismatches |
| Thm 8.12 | `quadric/dim4_m2_cases.py 3 all`, `5 all s 4` (4 shards) | 9,801 and 488,125 subspaces; 216 of type E2 (q = 9), 15,000 of type E1 (q = 25); no failure |
| Thm 8.12, colourings E1/E2 | `quadric/e1e2_slicelinear_check.py 3 all` | 216 checked, 216 verified |
| Thm 8.12, colourings E1/E2, q = 25 | `quadric/e1e2_slicelinear_check.py 5 all` | 15,000 checked, 15,000 verified |
| Cor 8.13, q = 9 | `quadric/dim_m2_check.py 3 4` | 11,011 subspaces; 648 with all weights 2; 0 counterexamples |
| Remark 8.14 (q = 27, periodic colourings) | `quadric/cross_quotient.py 3 3 1 1`, `3 3 2 1`, `3 3 2 2`; `quadric/periodic_search.py 3 3 2 3` | no tiling: 266 and 320 subspaces W of dimension 1 (a = 1, 2), 5,016 of dimension 2 (a = 2); 149,880 subspaces U of dimension 3 (a = 2) |
| trivial stabilizers | `stabilizer_check.py` | all trivial |
| Thm 9.2 (order 12) | `h12.py`, `h12b.py`, `h12c.py` | all 15,400 triple partitions give dim 8; dimension formula of (b) observed |
| Thm 9.6, Thm 9.8 | `h48_check.py 23`, `h48_check.py 47` | orbit partitions give dim K = 8 and 16; slices are the orbits |
| Lemma 9.7 | `case_check.py` | 1,112 profiles, all D ≤ 8 |
| Thm 9.8 consistency | `falsify.py 47 30` | no D above 8 |
| remark after Thm 9.8 | `sextet.py 23` | D = 5, dim K = 7 |
| §9, order-24 KS pair | `check24.py` | 256 orthogonal bases, 2,048 rays, each vector in 3 bases |
| Table 2 | `verify_table2_witnesses.py` | the stored witnesses for orders 12, 20, 24, 28, 36 reproduce the rows of the matrices, dim K = 8, 16, 8, 24, 24 and the slice weights; orders 44 and 48 follow from Thm 9.2(e) and Thm 9.8 |
| Prop 9.3 | `paley_rank_check.py` | all ranks match |
| Remark 11.7 (SAT) | `line_conj_sat.py 7 2 --prune` | 517 configurations pruned to 12 up to symmetry; all 12 UNSAT, none satisfiable (5,928 s) |
| Remark 11.7 (classification), line-slice evidence | `butson_check.py butson/BH-10-5.txt butson/BH-14-7.txt` | 12,348 directions per class of BH(14,7), each with a slice of size 1 or p |
| Prop 11.8, permutation planes | `rank_bound_check.py` | bound holds on 5,400 directions; 90 / 6 / 84, 7 planes of 18 |
| Prop 11.2 | `equidist_check.py` | 259,308 tests, no violation |
| Thm 12.1 (a) | `sloane16.py sloane/had.16.*.txt` | every normalization tiles |
| Thm 12.1 (b) | `sub16.py sloane/had.16.{1,2,3,4}.txt` | every partition colourable |
| Thm 12.1 (c) | `ghz.py 3 3 30` | all colourable; no failure |
| Example 12.3 | `nontile_example.py` | determines all 6 directions; does not tile |
| Lean | `lake env lean lean/Conic.lean` (in a Lean v4.34.1 project with Mathlib) | compiles, no `sorry` |

## Rows checked against their saved logs

These runs take hours or need a cluster; their logs are in the archive and agree with the manifest.

| Manifest row | Log | Content |
|---|---|---|
| Thm 9.6 (order-24 search, C program) | `logs/h24_codim4.log`, `logs/h24_codim3.log` | 866,251,507 subspaces of codimension 4, 1,771 with all weights ≥ 3, none with all weights divisible by 3; control: 50,955,971 of codimension 3, 1,700,160 with all weights divisible by 3 |
| Cor 8.13, q = 25 | `quadric/logs/m2_dim4_q25.log` | 508,431 subspaces; 15,000 with all weights 2; 0 counterexamples |
| Thm 8.12, q = 49, 121, 169 | `quadric/logs/e1e2_p7.log`, `e1e2_p11.log`, `e1e2_p13.log` | 692, 150, 60 colourings verified |
| sentence after Thm 8.12 (q = 49) | `quadric/logs/dim4_m2_e2_q49_sample.log` | 692 of 692 colourable (485 with f a square, 207 with f constant) |

## Rows not rerun

`linear_col.py` (not cited in the paper), the line-slice evidence after Remark 11.7 (`p11_ghsearch.py`,
`quadric/dim1_all.py`; the remark reports what these searches found, and their logs are included), and
`quadric/minbases_check.py` for q = 27 (not cited in the paper).
