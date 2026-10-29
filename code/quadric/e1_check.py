"""Quadric family: check the classification of subspaces K containing e1 = (1,0,0,0).
Claim: B(Tr D_q, K) is a KS pair iff K is not in {z3=0} and K0 = K cap {z3=0} projects onto F_q in z2 and in z4;
then every slice has p^k rows P(x) and p^k rows R(y), k = dim K - m - 2, and the counting equation has no solution.
Otherwise B is colourable (checked here by exact cover when no slice has weight 1).
Rows in the coordinates of Thm 8.2: P(x) = (Q1(x), x, 0, 0), R(y) = (Q2(y), 0, 1, y), Q1 = Tr(cx^2), Q2 = Tr(cx^2/n), c = (n-1)/4."""
import sys, io, contextlib, time
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
    from quadric_gh import field
from quadric_ks import dual_subspaces, count_ok, colorable, algx
import itertools
def single_tile(A, slice_rows):
    """exact cover of K by translates of one tile (differences of the rows of one slice)"""
    K = [v for v in itertools.product(range(p), repeat=r) if all(sum(a * b for a, b in zip(f, v)) % p == 0 for f in A)]
    idx = {k: i for i, k in enumerate(K)}
    T = [tuple((a - b) % p for a, b in zip(slice_rows[0], h)) for h in slice_rows]
    Y = {t: [idx[tuple((a + b) % p for a, b in zip(t, s))] for s in T] for t in K}
    X = {j: set() for j in range(len(K))}
    for key, cols in Y.items():
        for j in cols: X[j].add(key)
    return algx(X, Y, [])
p, m = int(sys.argv[1]), int(sys.argv[2])
COLMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 10**9      # exact cover only when |K| <= COLMAX
F = field(p, m); q = p ** m; r = 2 * m + 2
nsq = next(a for a in F.els if a != F.zero and not F.is_square(a)); ninv = F.inv(nsq)
c = F.mul(F.sub(nsq, F.one), F.inv(F.const(4)))          # c = (n-1)/4 as in quadric_gh.build
Q1 = lambda x: F.trace(F.mul(c, F.mul(x, x))); Q2 = lambda y: F.trace(F.mul(F.mul(c, ninv), F.mul(y, y)))
Z = tuple([0] * m)
rows = [(Q1(x),) + x + (0,) + Z for x in F.els] + [(Q2(y),) + Z + (1,) + y for y in F.els]
kind = ['P'] * q + ['R'] * q

def rank(M):
    M = [list(v) for v in M]; rk = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]; inv = pow(M[rk][c], p - 2, p)
        M[rk] = [v * inv % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
def nullspace(A):          # basis of {v : A v = 0} in F_p^r
    M = [list(f) for f in A]; piv = []; rk = 0
    for c in range(r):
        pr = next((i for i in range(rk, len(M)) if M[i][c] % p), None)
        if pr is None: continue
        M[rk], M[pr] = M[pr], M[rk]; inv = pow(M[rk][c], p - 2, p); M[rk] = [v * inv % p for v in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]; M[i] = [(a - f * b) % p for a, b in zip(M[i], M[rk])]
        piv.append(c); rk += 1
    basis = []
    for fc in [c for c in range(r) if c not in piv]:
        v = [0] * r; v[fc] = 1
        for i, c in enumerate(piv): v[c] = (-M[i][fc]) % p
        basis.append(v)
    return basis

z3 = [0] * r; z3[m + 1] = 1
stats = Counter(); bad = 0; t0 = time.time()
for k in range(0, r):                                   # annihilator dimension; f(e1) = 0 for all f
    for A0 in dual_subspaces(p, r - 1, k):
        A = [[0] + f for f in A0]; dimK = r - k
        mixed_possible = rank(A + [z3]) > rank(A) if k else True     # K not inside {z3 = 0}
        K0 = nullspace(A + [z3]) if k else nullspace([z3])
        u2 = rank([v[1:m + 1] for v in K0]) if K0 else 0
        u4 = rank([v[m + 2:] for v in K0]) if K0 else 0
        predicted_ks = mixed_possible and u2 == m and u4 == m
        lab = {}; members = {}
        for z, kd in zip(rows, kind):
            key = tuple(sum(a * b for a, b in zip(f, z)) % p for f in A)
            lab.setdefault(key, Counter())[kd] += 1
            members.setdefault(key, []).append(z)
        ws = sorted(sum(c.values()) for c in lab.values())
        ks = not count_ok(ws, p ** dimK)
        if predicted_ks:
            kk = dimK - m - 2
            ok = ks and kk >= 0 and all(c['P'] == p ** kk and c['R'] == p ** kk for c in lab.values())
            stats[('KS', dimK)] += 1
            if not ok: bad += 1; print('MISMATCH (predicted KS):', A, ws)
        else:
            if ks: bad += 1; print('MISMATCH (counting fails but predicted colourable):', A, ws); continue
            if 1 in ws: stats[('colourable: weight-1 slice', dimK)] += 1; continue
            pure = [key for key, c in lab.items() if len(c) == 1]        # slices with rows of one kind only
            if any(single_tile(A, members[key]) for key in pure[:3]):
                stats[('colourable by one pure tile', dimK)] += 1; continue
            if p ** dimK > COLMAX: stats[('no weight 1, not tested', dimK)] += 1; continue
            col = colorable(rows, A, p, r)
            stats[('colourable by exact cover' if col else 'NOT COLOURABLE', dimK)] += 1
            if not col: bad += 1; print('NOT COLOURABLE:', A, ws)
print(f'q={p}^{m}: ' + '; '.join(f'{k[0]} dim {k[1]}: {v}' for k, v in sorted(stats.items(), key=lambda t: (t[0][1], t[0][0]))))
print('mismatches:', bad, f'[{time.time() - t0:.0f}s]')
