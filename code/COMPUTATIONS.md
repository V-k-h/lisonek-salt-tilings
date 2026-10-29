# Computations: manifest

Each computational statement in `paper.tex`, with the script that produces it, its inputs, the number of cases checked, and the outcome.
Scripts are run from `code/` (or `code/quadric/` where stated) with Python 3. The SAT scripts use PySAT (CaDiCaL 1.5.3).

| Statement in the paper | Script | Inputs, normalizations | Cases | Outcome |
|---|---|---|---|---|
| Thm 4.1: E-basis identities and that D is a GH(Z_p,2) | `p3mod4_proof_check.py` | Construction 1, primes p ≡ 3 (mod 4) with p ≤ 47 | 8 primes | all identities hold |
| Cor 5.2: 15, 45, 91 KS hyperplanes for p = 3, 5, 7 | `verify.py` | Jungnickel D | all hyperplanes, p ≤ 13 | counts p(2p−1) |
| Thm 6.1: dim K ≤ 2 has a weight-1 slice | `verify.py` | Jungnickel D, p ≤ 13 | all subspaces of dim 1, 2 | none without a weight-1 slice |
| After Thm 7.5: colourings of Thm 7.4 tested against every basis | `p3mod4_proof_check.py` | class 2b(iii) hyperplanes | p = 3, 7, 11: every basis (cited); p = 19, 23: 15 sampled hyperplanes each; all p ≡ 3 mod 4 up to 71 through the values of β (not cited) | all valid |
| (not cited in the paper) colourings with \|C_v\| ≤ 1 | `linear_col.py` (logs in `logs/`) | p = 7, 11, 19 | 148, 606, 3250 hyperplanes | all found and checked against every basis |
| Prop 8.1, Thm 8.2: quadric family is GH, rank 2m+2 | `quadric/quadric_gh.py` | q ∈ {3, 5, 9, 25, 27, 49} | — | GH, rank 2m+2 |
| Thm 8.5: hyperplane weight classes and counts | `quadric/quadric_weights.py` | q ∈ {3, 9, 25, 27} | all hyperplanes | as stated |
| Thm 8.10 (m = 2): one-plane colourings | `quadric/separable.py` | q ∈ {9, 25} | all hyperplanes | all colourable |
| Remark after Prop 8.9: dim K ≤ m+1 colourable (Thm 8.8), and the subspaces of Prop 8.9 (both proved; checked as well; Prop 8.9 also by `quadric/e1_check.py`: 16 and 36 in dim 3) | `quadric/minbases_check.py` | q = 9, 25 | 45,255 and 3,070,893 subspaces | all colourable; 16 and 36 without a weight-1 slice |
| (not cited in the paper) Thm 8.11 (proved): every subspace containing e1; KS iff all rows mixed, then slices p^k P-rows + p^k R-rows; otherwise one pure slice tiles K; q·|GL_m(F_p)| KS subspaces of dim m+2 | `quadric/e1_check.py p m 0` (log `quadric/logs/e1_check.log`) | rows in the coordinates of Thm 8.2 with c = (n−1)/4; q ∈ {3, 5, 7, 9, 25, 27} | all subspaces containing e1 (2664 for q = 9, 42,176 for q = 25, 2,052,656 for q = 27) | 0 mismatches; KS of dim m+2: 6, 20, 42, 432, 12000, 303264 (= q·|GL_m(F_p)|); every non-KS subspace without a weight-1 slice is tiled by one pure tile (exact cover) |
| (not cited in the paper) Thm 8.12 (proved): the case analysis for m = 2, dim K = 4, e1 ∉ K. A1, A2 and f are computed from the formulas of the proof; checks χ(A1A2) = −ε1ε2, the vertex values against the slices, every predicted weight-1 / all-weight-2 / exceptional class, and that all weights are 2 only for opposite signs with f = 0 | `quadric/dim4_m2_cases.py p all` (p = 5 in 4 shards: `5 all s 4`), `p sample n`; `--sat N` checks colourability | p = 3, 5 all; p = 7 (30,000), 11 and 13 (20,000 each) sampled | 9801 and 488,125 subspaces with e1 ∉ K | no failure; all-weight-2: 216 and 3000 (= p³(p²−1)); exceptional: 216 of type E2 with f constant (q = 9, all colourable by SAT) and 15,000 of type E1 (q = 25) |
| (consistency check of Thm 8.12, which is proved) the colourings of types E1, E2 in its proof | `quadric/e1e2_slicelinear_check.py p all`, `... p sample n max` (logs `quadric/logs/e1e2_p{3,5,7,11,13}.log`) | m = 2, exceptional 4-dim subspaces (types E1, E2) as classified by `dim4_m2_cases.classify`; q = 9, 25: all subspaces; q = 49: 30,000 random; q = 121, 169: 40,000 random, first 150 resp. 60 exceptional ones checked | 216 (E2), 15,000 (E1), 692 (E2), 150 (E2), 60 (E1) | every structural claim of the proof holds; every colouring meets each basis exactly once |
| Sentence after Thm 8.12: type E2 for q = 49 (sample) | `quadric/dim4_m2_e2_sample.py` (log `quadric/logs/dim4_m2_e2_q49_sample.log`) | 30,000 random 4-dim subspaces without e1, seed 7 | 692 of type E2 (485 with f a square, 207 with f constant) | all colourable (SAT, models checked); f constant takes 1–10 min each, f square seconds |
| Cor 8.13 and the paragraph after Thm 8.12: every 4-dim subspace for q = 9 and q = 25 | `quadric/dim_m2_check.py 3 4` (q = 9, Algorithm X); `quadric/dimcheck_cluster.py dimcheck 5 4 --nshards 100` then `merge 5 4` (q = 25, exact cover by CaDiCaL with each model checked; log `quadric/logs/m2_dim4_q25.log`) | rows in the coordinates of Thm 8.2 | 11,011 and 508,431 subspaces | 0 counterexamples; KS by counting 648 and 15,000 (all weights 2); q = 25 exceptional slice weights (2^9, 3^8, 4^2) |
| (not cited in the paper) q = 27, dim ≤ 3 | `quadric/minbases_check.py` | q = 27 | 25,994,820 subspaces | all colourable; dim 4 not completed |
| (consistency check; Lemma 8.3, Thm 9.2 and the remark after Thm 9.8 prove these) trivial translation stabilizers | `stabilizer_check.py` | D (p ≤ 11), quadric q ∈ {9, 25, 27}, normalized Paley orders 12 and 24 | — | all trivial |
| (consistency check of Thm 9.2, which is proved) order 12: rank 10; all 15,400 triple partitions give dim 8 | `h12.py`, `h12b.py`, `h12c.py` | normalized Paley order 12 | 15,400 partitions | dim 8 |
| Thm 9.6 (proved) and the remark after it: the 1771 dimension-7 subspaces with all slice weights ≥ 3 | `h24_min.py`, `h24_min.c` (`h24_min 4 < h24_rows.txt`; control `h24_min 3`) | normalized Paley order 24, all subspaces of row(R) of codimension 4 (resp. 3) | 866,251,507 (resp. 50,955,971) | 1771, all with weights 4^6, none with all weights divisible by 3 (control: 1,700,160 with all weights divisible by 3) |
| Thm 9.6: the orbit partition of x ↦ −1/(x+1) gives dim K = 8 (proved via Lemma 9.4); the slices are the eight orbits (now proved via the self-duality in Lemma 9.4; checked as well); Col(R)+⟨1⟩ has the Golay weight distribution (consistency check) | inline check, see `logs/h24_codim4.log`; again `h48_check.py 23` (log `logs/h48_check.log`) | normalized Paley order 24 | — | dim 8, slices = orbits; weights 1, 759, 2576, 759, 1 |
| Thm 9.8: G = Col(R)+⟨1⟩ is a doubly even self-dual [48,24,12] code (consistency check of the cited fact); the 16 orbits of x ↦ −1/(x+1) give dim K = 16, and the slices are the orbits (now proved; checked as well) | `h48_check.py 47` (log `logs/h48_check.log`) | normalized Paley order 48 | all 2^24 codewords | weights 1, 17296, 535095, 3995376, 7681680, 3995376, 535095, 17296, 1; dim K = 16; slices = orbits |
| (not cited in the paper) Lemma 9.7, the case analysis in its proof: floor(o/2) + bound on dim C_E ≤ 8 | `case_check.py` (log `logs/h48_cases.log`) | every multiset of even part sizes with the worst odd count | 1112 profiles | all ≤ 8 with bounds (i) and (ii); no exact search needed |
| (not cited in the paper) Thm 9.8 consistency: D for orbit partitions of fixed-point-free elements of PSL_2(47), plus a local search | `falsify.py 47 30` (log `logs/h48_falsify.log`) | normalized Paley order 48 | cycle types 3^16, 4^12, 6^8, 8^6, 12^4, 24^2; 30 local searches | D = 8, 6, 4, 3, 2, 1; none above 8 |
| (not cited in the paper) remark after Thm 9.8: a Golay sextet has D = 5, dim K = 7 | `sextet.py 23` (log `logs/h24_sextet.log`) | normalized Paley order 24 | one sextet | D = 5 |
| §9: order-24 KS pair checked directly | `check24.py` | normalized Paley order 24 | 256 bases | orthogonal, 2048 rays, incidence 3 |
| Table 2, order 24 (orders 12, 20, 28, 36, 44 follow from Thm 9.2(e); order 48 is Thm 9.8, the earlier randomized value 2^20 is superseded) | `modq.py`, `table1_witnesses.py` (search; it enumerates K, so orders 44 and 48 need more memory than a desktop); `verify_table2_witnesses.py` (checks the stored witnesses) | Paley I q ∈ {11, 19, 23, 43, 47}, Paley II q ∈ {13, 17}, normalized | randomized search | witnesses in `witnesses/table1_partitions.json` for orders 12–36; order 48 witness is the orbit partition of Thm 9.8 |
| (consistency check) Prop 9.3 and the normalized ranks after it | `paley_rank_check.py` | Paley I: 13 primes; Paley II: 8 primes | — | all match; matrix identity holds |
| Remark 11.7: λ = 2, p = 7 (SAT) | `line_conj_sat.py 7 2 --prune` (log `lineconj_p7_l2_pruned.log`) | 12 gap classes, zero row in every slice and position | 12 classes | all UNSAT. DRAT certificates: generator `make_certs.py`, not yet produced |
| Remark 11.7: λ = 2, p = 7 (classification) | `butson_check.py` | BH(14,7), 3 classes, all 196 normalizations | 12,348 directions each | no all-weight-2 direction |
| Line-slice evidence (after Remark 11.7) | `p11_ghsearch.py`, `quadric/dim1_all.py`, `butson_check.py` | Jungnickel p = 11 and transpose (all normalizations), Kronecker sums (5 normalizations each), quadric q ∈ {25, 27, 49, 121}, BH(10,5), BH(14,7) | see logs | no exception |
| (not cited in the paper) Prop 11.8 checked | `rank_bound_check.py` | BH(14,7) (normalizations with i0, j0 ≡ 0 mod 3), quadric q = 25 | 5,400 directions | bound holds |
| p = 7 permutation planes (after Prop 11.8) | `rank_bound_check.py` | — | 720 permutations | 90 / 6 / 84, 7 planes of 18 |
| (not cited in the paper) Prop 11.2 checked | `equidist_check.py` | BH(14,7), all normalizations | 259,308 tests | no violation |
| Thm 12.1 (a) | `sloane16.py` | 5 Sloane classes of order 16, all normalizations | 16 + 64 + 128 + 256 + 128 row sets | all tile |
| Thm 12.1 (b) | `sub16.py` | 4 non-Sylvester classes, first 4 normalized row sets in lexicographic order | all subspace partitions | all colourable |
| Thm 12.1 (c) | `ghz.py 3 3 30` | 30 random GH(Z_3,3), `random.Random(7)` (fixed seed) | all submodules spanned by row differences | all colourable by exact cover |
| Example 12.3 (non-tiling 5-set) | `nontile_example.py` | — | all 5-subsets as T | does not tile |
| (not cited in the paper) Lean: n⁻¹w₄² − w₂² ≠ 0 for n a nonsquare, used in Thm 6.1 | `lean/Conic.lean` | Mathlib v4.34.1 | — | compiles without `sorry` |

Data: Sloane order-16 matrices in `sloane/`, Butson matrices in `butson/`.

## Class C hyperplanes for q = 27 (Remark 8.14; open)

| Question | Script | Outcome |
|---|---|---|
| reduction: every class C hyperplane is equivalent to ker(z1 + a z3); a = 1 weights (15,15,24), a = 2 weights (21,12,21) | `quadric/cross_tiling.py` (the reduction is proved in Remark 8.14) | two hyperplanes; O(Q1), O(Q2) of order 48 |
| cross tilings of F_27^2 periodic under W | `quadric/cross_quotient.py` | none: dim W = 1 (a = 1, 2), dim W = 2 (a = 2) |
| colourings periodic under U ≤ K | `quadric/periodic_search.py` | none for a = 2, dim U = 3 (149,880 U); dim U = 2 timed out at 50 min with none found |
