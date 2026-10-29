"""Quadric family, m = 2 (q = p^2), dim K = 4, e1 not in K: check the case analysis of the m = 2 dimension-4 theorem.
K = ker mu cap ker nu with mu(e1) = 0, nu(e1) = 1, mu = (0, a2, a3, a4), nu = (1, b2, b3, b4), lambda_1 = <a2,.>,
lambda_2 = <a4,.> (coordinates z1 | z2 | z3 | z4 of Theorem 8.2).  Cases:
  1   K in {z3=0} (a2 = a4 = 0)                         -> colourable (level sets)
  2A  exactly one of a2, a4 is 0                          -> colourable (level sets / subgroup tiles)
  2B  a2, a4 != 0, W_i = ker lambda_i, eps_i = chi(Q_i(w_i)):
      isotropic W_1 or W_2                                -> weight-1 slice
      non-isotropic: A1 = Q1(x1), A2 = Q2(y2) (x1 in W1^perp, lambda_1(x1) = 1, etc.), f = m1 - m2 computed from the
      formulas of the proof; predictions: weight-1 slice / all weights 2 (iff opposite signs and f = 0) / exceptional E1, E2.
For every K the predicted vertex values are compared with the slices, and the prediction with the slice weights;
--sat N checks colourability by SAT for up to N subspaces per case (cases 1, 2A, E1, E2).
Usage: dim4_m2_cases.py p (all [shard nshards] | sample n) [--sat N]"""
import sys, itertools, random, time
from collections import Counter
from cross_tiling import forms
import quadric_gh

quadric_gh.FIELDS.update({(13, 2): [2, 0], (17, 2): [3, 0]})        # t^2 = 2, 3: non-squares mod 13, 17

def leg(a, p): a %= p; return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)

def colorable_sat(rows, A, p):
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Cadical153
    K = [v for v in itertools.product(range(p), repeat=6) if all(sum(a * b for a, b in zip(f, v)) % p == 0 for f in A)]
    idx = {k: i for i, k in enumerate(K)}; groups = {}
    for z in rows: groups.setdefault(tuple(sum(a * b for a, b in zip(f, z)) % p for f in A), []).append(z)
    pls = []; cover = [[] for _ in K]
    for g, rs in groups.items():
        T = [tuple((a - b) % p for a, b in zip(rs[0], h)) for h in rs]
        for t in K:
            cells = [idx[tuple((a + b) % p for a, b in zip(t, s))] for s in T]; pls.append(cells)
            for c in cells: cover[c].append(len(pls))
    top = len(pls); cls = []
    for lits in cover:
        enc = CardEnc.equals(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter); top = max(top, enc.nv); cls += enc.clauses
    with Cadical153(bootstrap_with=cls) as S:
        if not S.solve(): return False
        mod = S.get_model()
    cnt = [0] * len(K)
    for i, cells in enumerate(pls):
        if mod[i] > 0:
            for c in cells: cnt[c] += 1
    assert all(c == 1 for c in cnt)
    return True

class Setup:
    def __init__(self, p):
        self.p = p; F, Q1, Q2 = forms(p, 2); self.F, self.Q1, self.Q2 = F, Q1, Q2
        self.els = F.els
        self.rows = [(Q1[x],) + x + (0,) + (0, 0) for x in F.els] + [(Q2[y],) + (0, 0) + (1,) + y for y in F.els]
        self.i2, self.i4 = pow(2, -1, p), pow(4, -1, p)
    def B(self, Q, x, y):
        p = self.p; return (Q[self.F.add(x, y)] - Q[x] - Q[y]) * self.i2 % p

def dot(a, x, p): return (a[0] * x[0] + a[1] * x[1]) % p

def classify(S, nu, mu):
    """case label, predicted property ('w1', 'all2', 'tile', 'E1', 'E2', 'p3'), and predicted (m1, m2) or None"""
    p = S.p; a2, a3, a4 = mu[1:3], mu[3], mu[4:6]; b2, b3, b4 = nu[1:3], nu[3], nu[4:6]
    z2, z4 = not any(a2), not any(a4)
    if z2 and z4: return "1", "tile", None
    if z2 or z4: return "2A", "tile", None
    w1 = (a2[1] % p, -a2[0] % p); w2 = (a4[1] % p, -a4[0] % p)          # kernels of lambda_1, lambda_2
    e1, e2 = leg(S.Q1[w1], p), leg(S.Q2[w2], p)
    if e1 == 0 or e2 == 0:
        assert not (e1 == 0 and e2 == 0)                                  # one of Q1, Q2 is anisotropic
        return "2B isotropic", "w1", None
    x1 = next(x for x in S.els if S.B(S.Q1, x, w1) == 0 and dot(a2, x, p)); x1 = tuple(c * pow(dot(a2, x1, p), -1, p) % p for c in x1)
    y2 = next(y for y in S.els if S.B(S.Q2, y, w2) == 0 and dot(a4, y, p)); y2 = tuple(c * pow(dot(a4, y2, p), -1, p) % p for c in y2)
    A1, A2 = S.Q1[x1], S.Q2[y2]
    assert leg(A1 * A2, p) == -e1 * e2                                    # opposite discriminants (Lemma 8.4)
    m1 = lambda u: (A1 * u * u + dot(b2, x1, p) * u - dot(b2, w1, p) ** 2 * S.i4 * pow(S.Q1[w1], -1, p)) % p
    m2 = lambda u: (A2 * (u - a3) ** 2 + dot(b4, y2, p) * (u - a3) + b3 - dot(b4, w2, p) ** 2 * S.i4 * pow(S.Q2[w2], -1, p)) % p
    fv = [(m1(u) - m2(u)) % p for u in range(p)]
    c0 = fv[0]; a = (fv[2] - 2 * fv[1] + c0) * S.i2 % p; b = (fv[1] - c0 - a) % p          # f = a u^2 + b u + c0
    assert all((a * u * u + b * u + c0 - fv[u]) % p == 0 for u in range(p))
    assert a == (A1 - A2) % p
    eps = e1; vert = ([m1(u) for u in range(p)], [m2(u) for u in range(p)])
    if e1 == e2:
        lab = "2B same"
        if p % 4 == 3: return lab + " p=3(4)", "w1", vert
        if (b * b - 4 * a * c0) % p: return lab + " p=1(4) Delta!=0", "w1", vert
        return (lab + " p=1(4) Delta=0 -", "w1", vert) if eps * leg(a, p) == -1 else ("E1", "E1", vert)
    lab = "2B opposite"
    if not any(fv): return lab + " f=0", "all2", vert
    if p % 4 == 1: return lab + " p=1(4) f!=0", "w1", vert
    if a == 0 and b: return lab + " p=3(4) deg f=1", "w1", vert
    if a == 0: return (lab + " p=3(4) f=c, eps chi(c)=+1", "w1", vert) if eps * leg(c0, p) == 1 else ("E2 (f constant)", "E2", vert)
    if (b * b - 4 * a * c0) % p: return (lab + " p=3(4) Delta!=0", "w1", vert) if p >= 7 else ("p=3, Delta!=0", "p3", vert)
    return (lab + " p=3(4) Delta=0 +", "w1", vert) if eps * leg(a, p) == 1 else ("E2 (f square)", "E2", vert)

def check(S, A, stats, bad, satn, satc):
    p = S.p; r0, r1 = A
    if r0[0] % p == 0: r0, r1 = r1, r0
    if r0[0] % p == 0: return                                             # e1 in K
    nu = [x * pow(r0[0], -1, p) % p for x in r0]; mu = [(x - r1[0] * y) % p for x, y in zip(r1, nu)]
    case, pred, vert = classify(S, nu, mu)
    lab = [(sum(a * b for a, b in zip(mu, r)) % p, sum(a * b for a, b in zip(nu, r)) % p) for r in S.rows]
    q = p * p; c1 = Counter(lab[:q]); c2 = Counter(lab[q:]); wt = c1 + c2
    has1 = 1 in wt.values(); all2 = set(wt.values()) == {2}
    if vert is not None:                                                  # vertex values of the proof match the slices
        for u in range(p):
            if c1[(u, vert[0][u])] != 1 or c2[(u, vert[1][u])] != 1: bad["vertex values differ from the slices"] += 1; break
    stats[case] += 1
    ok = {"w1": has1, "all2": all2, "tile": not all2, "E1": not has1 and not all2, "E2": not has1 and not all2, "p3": True}[pred]
    if not ok: bad[f"{case}: prediction {pred} fails"] += 1
    if all2 and pred != "all2": bad["all weights 2 outside 'opposite, f = 0'"] += 1
    if pred in ("tile", "E1", "E2", "p3") and not has1 and satc[case] < satn:
        satc[case] += 1; stats[f"{case}: colourable by SAT = {colorable_sat(S.rows, [nu, mu], p)}"] += 1

def main():
    p = int(sys.argv[1]); mode = sys.argv[2]
    satn = int(sys.argv[sys.argv.index("--sat") + 1]) if "--sat" in sys.argv else 0
    S = Setup(p); stats, bad, satc = Counter(), Counter(), Counter(); t0 = time.time(); tot = 0
    if mode == "all":
        shard, nsh = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 and sys.argv[3] != "--sat" else (0, 1)
        idx = -1
        for j in range(1, 6):                                             # RREF annihilators with pivots (0, j)
            free = [(0, c) for c in range(1, 6) if c != j] + [(1, c) for c in range(j + 1, 6)]
            for vals in itertools.product(range(p), repeat=len(free)):
                idx += 1
                if idx % nsh != shard: continue
                A = [[0] * 6, [0] * 6]; A[0][0] = 1; A[1][j] = 1
                for (i, c), v in zip(free, vals): A[i][c] = v
                tot += 1; check(S, A, stats, bad, satn, satc)
    else:
        rnd = random.Random(p); n = int(sys.argv[3])
        while tot < n:
            A = [[1] + [rnd.randrange(p) for _ in range(5)], [0] + [rnd.randrange(p) for _ in range(5)]]
            if not any(A[1]): continue
            tot += 1; check(S, A, stats, bad, satn, satc)
    print(f"p={p} {mode} {' '.join(sys.argv[3:])}: {tot} subspaces with e1 not in K [{time.time()-t0:.0f}s]")
    for k, v in sorted(stats.items()): print(f"   {k}: {v}")
    print("   FAILED:", dict(bad) if bad else "none", flush=True)

if __name__ == "__main__":
    main()
