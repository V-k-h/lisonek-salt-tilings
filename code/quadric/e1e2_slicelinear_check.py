"""Check the slice-linear colouring proposed for the exceptional subspaces (types E1, E2) of the m = 2 dimension-4
theorem.  For each exceptional K (classified by dim4_m2_cases.classify) this follows the proposed proof:
  * in a fibre mu = u with f(u) != 0, the vertex slices (u, m1(u)), (u, m2(u)) have weight 3 and the forms
    {P, R, R} and {P, P, R};
  * the fibre has (p-1)/2 (E1) or (p-3)/2 (E2) slices of weight 2, each with both rows in the same block;
  * with C0 such a slice and C3 the vertex slice whose pair lies in the other block, the tile differences
    d, e (of C3) and d0 (of C0) are linearly independent, and beta(d) = beta(e) = beta(d0) = 1 defines
    beta: K -> F_p (for p = 3: beta(d) = beta(e) = 1, no C0);
  * the translates w + A3 and w + t g + A0 (w in ker beta, t = 3, 5, ..., p - 2) partition K;
and then checks from the definition that Z = {r_C + translate} meets every basis R + k (k in K) exactly once.
Usage: e1e2_slicelinear_check.py p (all | sample n) [max_checked]"""
import sys, itertools, random, time
from collections import Counter
from dim4_m2_cases import Setup, classify

def span(basis, p):
    return {tuple(sum(c * b[i] for c, b in zip(cs, basis)) % p for i in range(6)): cs
            for cs in itertools.product(range(p), repeat=len(basis))}

def kernel(rows, p):
    """basis of {z in F_p^6 : r.z = 0 for r in rows}"""
    M = [list(r) for r in rows]; piv = []; i = 0
    for c in range(6):
        k = next((j for j in range(i, len(M)) if M[j][c] % p), None)
        if k is None: continue
        M[i], M[k] = M[k], M[i]; inv = pow(M[i][c], -1, p); M[i] = [x * inv % p for x in M[i]]
        for j in range(len(M)):
            if j != i and M[j][c] % p: M[j] = [(x - M[j][c] * y) % p for x, y in zip(M[j], M[i])]
        piv.append(c); i += 1
    basis = []
    for fc in (c for c in range(6) if c not in piv):
        v = [0] * 6; v[fc] = 1
        for r, pc in enumerate(piv): v[pc] = -M[r][fc] % p
        basis.append(tuple(v))
    return basis

def independent(vs, p):
    return not any(any(cs) and all(sum(c * v[i] for c, v in zip(cs, vs)) % p == 0 for i in range(len(vs[0])))
                   for cs in itertools.product(range(p), repeat=len(vs)))

def sub(x, y, p): return tuple((a - b) % p for a, b in zip(x, y))
def add(x, y, p): return tuple((a + b) % p for a, b in zip(x, y))

def check_one(S, nu, mu, typ, vert, bad):
    p, q = S.p, S.p * S.p; rows = S.rows; Rset = set(rows)
    sl = {}
    for i, r in enumerate(rows):
        key = (sum(a * b for a, b in zip(mu, r)) % p, sum(a * b for a, b in zip(nu, r)) % p)
        sl.setdefault(key, []).append((i < q, r))                        # (is P-row, row)
    m1, m2 = vert
    u = next((u for u in range(p) if (m1[u] - m2[u]) % p), None)
    if u is None: bad["no fibre with f(u) != 0"] += 1; return False
    s1, s2 = sl.get((u, m1[u]), []), sl.get((u, m2[u]), [])
    if sorted(isP for isP, _ in s1) != [False, False, True] or sorted(isP for isP, _ in s2) != [False, True, True]:
        bad["vertex slices not of the forms {P,R,R}, {P,P,R}"] += 1; return False
    two = [c for (uu, v), c in sl.items() if uu == u and len(c) == 2]
    if len(two) != ((p - 1) // 2 if typ == "E1" else (p - 3) // 2): bad["number of weight-2 slices in the fibre"] += 1; return False
    if any(c[0][0] != c[1][0] for c in two): bad["weight-2 slice with rows in both blocks"] += 1; return False
    basis = kernel([nu, mu], p); Kc = span(basis, p)
    if p == 3:
        C3 = s1; C0 = None
    else:
        C0 = two[0]; C3 = s1 if C0[0][0] else s2                           # C0 a P-pair -> use {P, R, R}
    single = next(r for isP, r in C3 if sum(1 for x, _ in C3 if x == isP) == 1)
    pair = [r for _, r in C3 if r != single]
    d, d2 = sub(single, pair[0], p), sub(single, pair[1], p); e = sub(d2, d, p)
    vecs = [d, e] + ([sub(C0[0][1], C0[1][1], p)] if C0 else [])
    if any(v not in Kc for v in vecs): bad["tile differences not in K"] += 1; return False
    co = [Kc[v] for v in vecs]
    if not independent(co, p): bad["d, e, d0 dependent"] += 1; return False
    beta = next(b for b in itertools.product(range(p), repeat=4) if all(sum(x * y for x, y in zip(b, c)) % p == 1 for c in co))
    bv = lambda k: sum(x * y for x, y in zip(beta, Kc[k])) % p
    W = [k for k in Kc if bv(k) == 0]; g = next(k for k in Kc if bv(k) == 1)
    tiles = [(single, w) for w in W]                                        # translate w of A3 <-> vector single + w
    if C0:
        r0 = C0[0][1]
        for t in range(3, p, 2):
            tg = tuple(t * x % p for x in g)
            tiles += [(r0, add(w, tg, p)) for w in W]
    A = {single: [sub(single, h, p) for _, h in C3]}
    if C0: A[r0] = [sub(r0, h, p) for _, h in C0]
    cov = Counter(add(tau, a, p) for rep, tau in tiles for a in A[rep])
    if len(cov) != len(Kc) or set(cov.values()) != {1} or any(k not in Kc for k in cov): bad["translates do not partition K"] += 1; return False
    Z = {add(rep, tau, p) for rep, tau in tiles}                           # the colouring, as vectors of R + K
    hits = Counter()
    for z in Z:
        for r in rows:
            k = sub(z, r, p)
            if k in Kc: hits[k] += 1
    if len(hits) != len(Kc) or set(hits.values()) != {1}: bad["Z is not a KS colouring"] += 1; return False
    return True

def main():
    p = int(sys.argv[1]); mode = sys.argv[2]
    cap = int(sys.argv[4]) if mode == "sample" and len(sys.argv) > 4 else (int(sys.argv[3]) if mode == "all" and len(sys.argv) > 3 else 10 ** 9)
    S = Setup(p); stats, bad = Counter(), Counter(); t0 = time.time(); tot = 0
    def gen():
        if mode == "all":
            for j in range(1, 6):
                free = [(0, c) for c in range(1, 6) if c != j] + [(1, c) for c in range(j + 1, 6)]
                for vals in itertools.product(range(p), repeat=len(free)):
                    A = [[0] * 6, [0] * 6]; A[0][0] = 1; A[1][j] = 1
                    for (i, c), v in zip(free, vals): A[i][c] = v
                    yield A
        else:
            rnd = random.Random(p); n = int(sys.argv[3]); k = 0
            while k < n:
                A = [[1] + [rnd.randrange(p) for _ in range(5)], [0] + [rnd.randrange(p) for _ in range(5)]]
                if any(A[1]): k += 1; yield A
    for r0, r1 in gen():
        tot += 1
        nu = [x * pow(r0[0], -1, p) % p for x in r0]; mu = [(x - r1[0] * y) % p for x, y in zip(r1, nu)]
        case, pred, vert = classify(S, nu, mu)
        if pred not in ("E1", "E2", "p3"): continue
        stats[f"{case}"] += 1
        if stats["checked"] >= cap: continue
        stats["checked"] += 1
        if check_one(S, nu, mu, "E1" if pred == "E1" else "E2", vert, bad): stats["colouring verified"] += 1
    print(f"p={p} {' '.join(sys.argv[2:])}: {tot} subspaces with e1 not in K [{time.time()-t0:.0f}s]")
    for k, v in sorted(stats.items()): print(f"   {k}: {v}")
    print("   FAILED:", dict(bad) if bad else "none", flush=True)

if __name__ == "__main__":
    main()
