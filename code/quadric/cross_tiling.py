"""Class C hyperplanes of the quadric family via cross tilings.
Up to the affine symmetries (x,y) -> (x+x0, y+y0), a class C hyperplane is K = ker(z1 + a z3), a != 0.
With beta(z) = (z2, z4), Proposition 7.1 (with beta to F_q^2) reduces colourability to tiling F_q^2 by
translates of the crosses T_v = (H_v x 0) u (0 x V_{v-a}), H_v = {Q1 = v}, V_w = {Q2 = w}; the colouring is
Z = {z : (z2, z4) in B_{alpha(z)}}.  This script finds such tilings by SAT and checks the colouring against
every basis R + k, k in K, in the coordinates of Theorem 8.2.
Usage: cross_tiling.py p m [a ...]"""
import sys, itertools, time
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
from quadric_gh import field

def forms(p, m):
    F = field(p, m)
    n = next(a for a in F.els if not F.is_square(a))
    c = F.mul(F.sub(n, F.one), F.inv(F.const(4))); ni = F.inv(n)
    sq = lambda x: F.mul(x, x)
    Q1 = {x: F.trace(F.mul(c, sq(x))) for x in F.els}
    Q2 = {y: F.trace(F.mul(F.mul(c, ni), sq(y))) for y in F.els}
    return F, Q1, Q2

def crosses(p, F, Q1, Q2, a):
    T = {}
    for v in range(p):
        H = [x for x in F.els if Q1[x] == v]
        V = [y for y in F.els if Q2[y] == (v - a) % p]
        cells = [(x, F.zero) for x in H] + [(F.zero, y) for y in V]
        assert len(set(cells)) == len(cells)          # the two arms are disjoint
        T[v] = cells
    return T

def solve(p, F, T, extra=None):
    G = [(s, u) for s in F.els for u in F.els]
    idx = {g: i for i, g in enumerate(G)}
    var, cover = {}, [[] for _ in G]
    for v, cells in T.items():
        for b in G:
            x = len(var) + 1; var[(v, b)] = x
            for (dx, dy) in cells:
                cover[idx[(F.add(b[0], dx), F.add(b[1], dy))]].append(x)
    top = len(var); cls = []
    for lits in cover:
        enc = CardEnc.equals(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter)
        top = max(top, enc.nv); cls += enc.clauses
    if extra: cls += extra(var)
    with Cadical153(bootstrap_with=cls) as S:
        if not S.solve(): return None
        mod = set(l for l in S.get_model() if l > 0)
    return {v: [b for b in G if var[(v, b)] in mod] for v in T}

def check_colouring(p, m, F, Q1, Q2, a, B):
    """Verify Z = {z : (z2,z4) in B_{alpha(z)}} meets every basis R+k, k in ker(z1 + a z3), exactly once."""
    rows = [(Q1[x], x, 0, F.zero) for x in F.els] + [(Q2[y], F.zero, 1, y) for y in F.els]
    Bs = {v: set(B[v]) for v in B}
    add = lambda z, k: ((z[0] + k[0]) % p, F.add(z[1], k[1]), (z[2] + k[2]) % p, F.add(z[3], k[3]))
    alpha = lambda z: (z[0] + a * z[2]) % p
    for z2 in F.els:
        for z4 in F.els:
            for z3 in range(p):
                k = ((-a * z3) % p, z2, z3, z4)
                hits = sum(1 for r in rows if (lambda z: (z[1], z[3]) in Bs[alpha(z)])(add(r, k)))
                if hits != 1: return False
    return True

if __name__ == "__main__":
    p, m = int(sys.argv[1]), int(sys.argv[2])
    F, Q1, Q2 = forms(p, m)
    for a in (map(int, sys.argv[3:]) if len(sys.argv) > 3 else range(1, p)):
        T = crosses(p, F, Q1, Q2, a)
        t0 = time.time(); B = solve(p, F, T)
        sizes = {v: len(T[v]) for v in T}
        if B is None:
            print(f"q={p**m} a={a}: cross sizes {sizes}: no tiling [{time.time()-t0:.0f}s]", flush=True); continue
        ok = check_colouring(p, m, F, Q1, Q2, a, B)
        print(f"q={p**m} a={a}: cross sizes {sizes}: tiling with {[len(B[v]) for v in sorted(B)]} translates; "
              f"colouring verified against all {p**(2*m+1)} bases: {ok} [{time.time()-t0:.0f}s]", flush=True)
