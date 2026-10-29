"""p=5, m=2 quadric family: for every hyperplane with solvable counting, test whether a
one-plane slice-linear colouring exists (a weight-p slice on which some functional beta is injective)."""
import itertools
from collections import Counter, defaultdict
from quadric_gh import build
from quadric_ks import dual_subspaces, count_ok
import sys; p,m=int(sys.argv[1]),int(sys.argv[2])
F,T,coords,_,_=build(p,m); r=len(coords[0])
duals=list(itertools.product(range(p),repeat=r))
dot=lambda a,z: sum(x*y for x,y in zip(a,z))%p
stats=Counter(); bad=[]
for A in dual_subspaces(p,r,1):
    a=A[0]; slices=defaultdict(list)
    for z in coords: slices[dot(a,z)].append(z)
    ws=tuple(sorted(len(v) for v in slices.values()))
    if not count_ok(ws,p**(r-1)): stats[("KS",ws)]+=1; continue
    found=False
    for v,pts in slices.items():
        if len(pts)!=p: continue
        diffs=[tuple((x-y)%p for x,y in zip(u,w)) for u,w in itertools.combinations(pts,2)]
        # beta must be non-parallel to alpha (beta(K)=F_p); injective on the slice
        for b in duals:
            if all(dot(b,d) for d in diffs):
                found=True; break
        if found: break
    if len(ws)==2 and ws[0]==ws[1]==p**m: found="transversal (Prop.)"
    stats[(ws,"one-plane colouring" if found else "NONE")]+=1
    if not found: bad.append(A)
print(dict(stats)); print("hyperplanes needing a different colouring:",len(bad))
