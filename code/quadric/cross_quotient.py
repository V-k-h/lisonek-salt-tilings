"""W-periodic cross tilings: tile G = F_q^2 by translates of the crosses T_v, with centre sets invariant under a
subspace W of G.  Such a tiling is the same as a tiling of G/W by the images of the crosses that are used, and
those images must be injective.  Enumerates all W of a given dimension, keeps the crosses that inject mod W,
and solves the quotient tiling by SAT.  Each tiling found is lifted to F_q^2 and the colouring is checked
against every basis with cross_tiling.check_colouring.
Usage: cross_quotient.py p m a dimW [max_solutions]"""
import sys, itertools, time
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
from cross_tiling import forms, crosses, check_colouring

def vec(F, g):            # (s, u) in F_q^2 -> tuple in F_p^{2m}
    return tuple(g[0]) + tuple(g[1])

def subspaces(p, n, k):
    """all k-dim subspaces of F_p^n as RREF bases"""
    for piv in itertools.combinations(range(n), k):
        free = [(i, j) for i in range(k) for j in range(n) if j > piv[i] and j not in piv]
        for vals in itertools.product(range(p), repeat=len(free)):
            B = [[0] * n for _ in range(k)]
            for i in range(k): B[i][piv[i]] = 1
            for (i, j), x in zip(free, vals): B[i][j] = x
            yield B

def quotient_map(p, n, B):
    """linear map F_p^n -> F_p^{n-k} with kernel span(B): coordinates at non-pivot columns after reduction"""
    k = len(B); piv = [next(j for j in range(n) if row[j]) for row in B]
    nonpiv = [j for j in range(n) if j not in piv]
    def q(v):
        v = list(v)
        for row, c in zip(B, piv):
            f = v[c]
            if f: v = [(x - f * y) % p for x, y in zip(v, row)]
        return tuple(v[j] for j in nonpiv)
    return q

def solve_quotient(p, cells, shapes):
    idx = {c: i for i, c in enumerate(cells)}
    var, cover = {}, [[] for _ in cells]
    add = lambda a, b: tuple((x + y) % p for x, y in zip(a, b))
    for v, S in shapes.items():
        for b in cells:
            x = len(var) + 1; var[(v, b)] = x
            for d in S: cover[idx[add(b, d)]].append(x)
    top = len(var); cls = []
    for lits in cover:
        enc = CardEnc.equals(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter)
        top = max(top, enc.nv); cls += enc.clauses
    with Cadical153(bootstrap_with=cls) as S:
        if not S.solve(): return None
        mod = set(l for l in S.get_model() if l > 0)
    return {v: [b for b in cells if var[(v, b)] in mod] for v in shapes}

if __name__ == "__main__":
    p, m, a, k = map(int, sys.argv[1:5]); maxsol = int(sys.argv[5]) if len(sys.argv) > 5 else 3
    F, Q1, Q2 = forms(p, m); T = crosses(p, F, Q1, Q2, a); n = 2 * m
    sizes = {v: len(T[v]) for v in T}
    print(f"q={p**m} a={a} cross sizes {sizes}; W of dimension {k}", flush=True)
    G = [(s, u) for s in F.els for u in F.els]
    cells = list(itertools.product(range(p), repeat=n - k))
    found = tried = 0; t0 = time.time()
    for B in subspaces(p, n, k):
        q = quotient_map(p, n, B)
        shapes = {}
        for v, S in T.items():
            img = [q(vec(F, c)) for c in S]
            if len(set(img)) == len(img): shapes[v] = img
        # counting: |G/W| must be a non-negative combination of the injective sizes
        N = p ** (n - k); szs = sorted({len(S) for S in shapes.values()})
        reach = {0}
        for t in range(1, N + 1):
            if any(t - s in reach for s in szs): reach.add(t)
        if N not in reach: continue
        tried += 1
        sol = solve_quotient(p, cells, shapes)
        if sol is None: continue
        # lift: centres in G whose image lies in the quotient centre set
        Bc = {v: [g for g in G if q(vec(F, g)) in set(sol.get(v, []))] for v in T}
        ok = check_colouring(p, m, F, Q1, Q2, a, Bc)
        found += 1
        print(f"  W = {B}: quotient tiling uses {[(v, len(sol[v])) for v in sorted(sol)]}; "
              f"lifted colouring verified: {ok} [{time.time()-t0:.0f}s]", flush=True)
        if found >= maxsol: break
    print(f"done: {tried} subspaces W passed the counting test, {found} tilings found [{time.time()-t0:.0f}s]", flush=True)
