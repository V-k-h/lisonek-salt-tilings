import io, contextlib, numpy as np
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, normalize_bits, rank
    from h24 import search
n=24; R=normalize_bits(paley1(23)); d,parts=search(R,(3,)*8,4000,30)
gens=[R[a]^R[P[0]] for P in parts for a in P[1:]]
basis=[]
for v in gens:
    for b in basis: v=min(v,v^b)
    if v: basis.append(v)
K={0}
for b in basis: K|={k^b for k in K}
print("dim K",len(basis),"|K|",len(K))
vec=lambda m: np.array([-1 if (m>>j)&1 else 1 for j in range(n)])
bases=[[r^k for r in R] for k in K]
ok=all(np.array_equal(np.array([vec(m) for m in B])@np.array([vec(m) for m in B]).T, n*np.eye(n,dtype=int)) for B in bases)
V=Counter(m for B in bases for m in B)
rays={min(m, m^((1<<n)-1)) for m in V}
print("all 256 bases orthogonal:",ok,"| vectors:",len(V),"| distinct rays:",len(rays),"| incidence per vector:",set(V.values()))
print("KS by counting: each vector in 3 bases, 3 ∤",len(K))
