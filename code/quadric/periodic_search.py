"""U-periodic colourings of a class C hyperplane K = ker(z1 + a z3) of the quadric family, for every subspace U of
K of dimension k.  In the coordinates (z2, z3, z4) of K, slice v has shape S_v = {(x,0,0): Q1(x)=v} u {(0,1,y):
Q2(y)=v-a}.  A U-periodic tiling of K is a tiling of K/U by the images of the shapes, which must be injective.
For each U passing the injectivity and counting tests the quotient tiling is solved by SAT; a solution gives the
colouring Z = {z : beta(z) in C_{alpha(z)}} with beta = quotient map on (z2, z3, z4), which is checked against
every basis R + k.
Usage: periodic_search.py p m a k [nproc]"""
import sys, itertools, time
import numpy as np
from multiprocessing import Pool
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
from cross_tiling import forms

P = M = A = Kd = None
SHAPES = None     # v -> int array (|S_v|, n) of shape vectors in F_p^n, n = 2m+1
ROWS = None

def setup(p, m, a):
    global P, M, A, SHAPES, ROWS, F, Q1, Q2
    P, M, A = p, m, a
    F, Q1, Q2 = forms(p, m)
    SHAPES = {}
    for v in range(p):
        S = [tuple(x) + (0,) + (0,) * m for x in F.els if Q1[x] == v] + \
            [(0,) * m + (1,) + tuple(y) for y in F.els if Q2[y] == (v - a) % p]
        SHAPES[v] = np.array(S, dtype=np.int64)
    # rows of R in V = F_p x F_q x F_p x F_q: (z1, z2, z3, z4); alpha = z1 + a z3
    ROWS = [(Q1[x], tuple(x), 0, (0,) * m) for x in F.els] + [(Q2[y], (0,) * m, 1, tuple(y)) for y in F.els]

def subspaces_with_pivot(n, k, piv, p):
    free = [(i, j) for i in range(k) for j in range(n) if j > piv[i] and j not in piv]
    for vals in itertools.product(range(p), repeat=len(free)):
        B = np.zeros((k, n), dtype=np.int64)
        for i in range(k): B[i, piv[i]] = 1
        for (i, j), x in zip(free, vals): B[i, j] = x
        yield B

def qmap(B, p):
    """matrix Qm (n x (n-k)) with v -> v @ Qm mod p a surjection with kernel rowspan(B) (B in RREF)"""
    k, n = B.shape
    piv = [int(np.nonzero(B[i])[0][0]) for i in range(k)]
    nonpiv = [j for j in range(n) if j not in piv]
    # v -> v - sum_i v[piv_i] B_i, then read nonpivot coordinates
    E = np.eye(n, dtype=np.int64)
    for i, c in enumerate(piv): E[c, :] = (E[c, :] - B[i, :]) % p
    return E[:, nonpiv] % p

def sat_tiling(p, s, shapes):
    cells = list(itertools.product(range(p), repeat=s))
    idx = {c: i for i, c in enumerate(cells)}
    var, cover = {}, [[] for _ in cells]
    for v, S in shapes.items():
        negS = [tuple((-x) % p for x in d) for d in S]
        for b in cells:
            x = len(var) + 1; var[(v, b)] = x
            for d in negS: cover[idx[tuple((y + z) % p for y, z in zip(b, d))]].append(x)
    top = len(var); cls = []
    for lits in cover:
        enc = CardEnc.equals(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter)
        top = max(top, enc.nv); cls += enc.clauses
    with Cadical153(bootstrap_with=cls) as S:
        if not S.solve(): return None
        mod = set(l for l in S.get_model() if l > 0)
    return {v: {b for b in cells if var[(v, b)] in mod} for v in shapes}

def check(Qm, C):
    p, m, a = P, M, A
    beta = lambda z2, z3, z4: tuple(int(x) for x in (np.array(z2 + (z3,) + z4) @ Qm) % p)
    for k2 in itertools.product(range(p), repeat=m):
        for k3 in range(p):
            for k4 in itertools.product(range(p), repeat=m):
                hits = 0
                for (r1, r2, r3, r4) in ROWS:
                    z2 = tuple((x + y) % p for x, y in zip(r2, k2)); z3 = (r3 + k3) % p
                    z4 = tuple((x + y) % p for x, y in zip(r4, k4))
                    z1 = (r1 - a * k3) % p
                    if beta(z2, z3, z4) in C.get((z1 + a * z3) % p, ()): hits += 1
                if hits != 1: return False
    return True

def work(args):
    p, m, a, k, piv = args
    setup(p, m, a); n = 2 * m + 1; s = n - k; N = p ** s
    sizes = sorted({len(S) for S in SHAPES.values()})
    out = []; tested = 0
    for B in subspaces_with_pivot(n, k, piv, p):
        Qm = qmap(B, p)
        shapes = {}
        for v, S in SHAPES.items():
            img = [tuple(r) for r in (S @ Qm) % p]
            if len(set(img)) == len(img): shapes[v] = img
        szs = sorted({len(S) for S in shapes.values()})
        reach = np.zeros(N + 1, bool); reach[0] = True
        for t in range(1, N + 1): reach[t] = any(t >= z and reach[t - z] for z in szs)
        if not reach[N]: continue
        tested += 1
        sol = sat_tiling(p, s, shapes)
        if sol is not None:
            out.append((B.tolist(), {v: len(c) for v, c in sol.items()}, check(Qm, sol)))
            break
    return tested, out

if __name__ == "__main__":
    p, m, a, k = map(int, sys.argv[1:5]); nproc = int(sys.argv[5]) if len(sys.argv) > 5 else 16
    setup(p, m, a); n = 2 * m + 1
    print(f"q={p**m} a={a} shape sizes {[len(SHAPES[v]) for v in range(p)]}; U of dimension {k} in K=F_{p}^{n}", flush=True)
    t0 = time.time(); tot = 0; found = []
    with Pool(nproc) as pool:
        for tested, out in pool.imap_unordered(work, [(p, m, a, k, piv) for piv in itertools.combinations(range(n), k)]):
            tot += tested
            for B, used, ok in out:
                found.append(B)
                print(f"  U = {B}: tiling uses {used}; colouring verified against all bases: {ok} [{time.time()-t0:.0f}s]", flush=True)
    print(f"done: {tot} subspaces U passed injectivity and counting; {len(found)} pivot classes with a tiling [{time.time()-t0:.0f}s]", flush=True)
