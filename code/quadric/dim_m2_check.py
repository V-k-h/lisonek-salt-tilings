"""Quadric family, m = 2 (q = p^2): test Conjecture 2.2 on every subspace K of dimension d of the row space,
in the coordinates of Theorem 8.2 (rows P(x) = (Q1(x), x, 0, 0) and R(y) = (Q2(y), 0, 1, y) in F_p^6).
For each K = ker A (A a (6-d) x 6 matrix in RREF): slice weights, counting condition, colourability by exact
cover when the counting condition holds and no slice has weight 1.  Also records the subspaces whose slices all
have weight 2, split by whether e1 = (1,0,0,0) lies in K, and by (k1, k2) and whether K lies in {z3 = 0}.
With shard i of n, only the subspaces with index = i (mod n) are checked, and a RESULT json line is printed
for merge_m2.py.
Usage: dim_m2_check.py p d [shard nshards]"""
import sys, itertools, time, json
from collections import Counter
import numpy as np
from cross_tiling import forms
from quadric_ks import count_ok, colorable

def main(p, d, shard=0, nsh=1):
    m = 2; r = 2 * m + 2; k = r - d
    F, Q1, Q2 = forms(p, m)
    rows = [(Q1[x],) + tuple(x) + (0,) + (0,) * m for x in F.els] + [(Q2[y],) + (0,) * m + (1,) + tuple(y) for y in F.els]
    Rm = np.array(rows, dtype=np.int64)
    stats = Counter(); all2 = Counter(); fails = []; t0 = time.time(); total = 0; idx = -1
    for pivs in itertools.combinations(range(r), k):
        free = [(i, j) for i in range(k) for j in range(r) if j > pivs[i] and j not in pivs]
        for vals in itertools.product(range(p), repeat=len(free)):
            idx += 1
            if idx % nsh != shard: continue
            A = np.zeros((k, r), dtype=np.int64)
            for i in range(k): A[i, pivs[i]] = 1
            for (i, j), v in zip(free, vals): A[i, j] = v
            total += 1
            lab = (Rm @ A.T) % p
            ws = sorted(Counter(map(tuple, lab)).values())
            e1in = not A[:, 0].any()
            z3sub = bool(np.linalg.matrix_rank(np.vstack([A, np.eye(r, dtype=np.int64)[m + 1]])) == k)  # K in {z3=0}?
            if set(ws) == {2}:
                all2[("e1 in K" if e1in else "e1 not in K", "K in z3=0" if z3sub else "K not in z3=0")] += 1
            if not count_ok(ws, p ** d):
                stats["KS by counting"] += 1; continue
            if 1 in ws:
                stats["solvable, weight-1 slice"] += 1; continue
            col = colorable([tuple(z) for z in rows], A.tolist(), p, r)
            stats[f"solvable, no weight 1, colourable={col}"] += 1
            if not col:
                fails.append((A.tolist(), ws)); print("   COUNTEREXAMPLE:", A.tolist(), ws, flush=True)
    tag = f"shard {shard}/{nsh} " if nsh > 1 else ""
    print(f"{tag}q={p}^2 dim K={d}: {total} subspaces [{time.time()-t0:.0f}s]")
    for key, n in sorted(stats.items()): print(f"   {key}: {n}")
    print("   all slice weights 2:", dict(all2))
    print("   counterexamples to Conjecture 2.2:", len(fails), flush=True)
    if nsh > 1:
        print("RESULT", json.dumps({"p": p, "d": d, "shard": shard, "nsh": nsh, "total": total, "secs": round(time.time() - t0),
                                   "stats": stats, "all2": {" | ".join(key): n for key, n in all2.items()}, "fails": fails}), flush=True)

if __name__ == "__main__":
    main(*map(int, sys.argv[1:5]))
