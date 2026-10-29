import sys
from collections import Counter
from quadric_gh import build
from quadric_ks import dual_subspaces, count_ok
for p,m in [(3,1),(3,2),(5,2),(3,3)]:
    F,T,coords,_,_=build(p,m); r=len(coords[0])
    hist=Counter()
    for A in dual_subspaces(p,r,1):
        f=A[0]; ws=tuple(sorted(Counter(sum(a*b for a,b in zip(f,z))%p for z in coords).values()))
        hist[ws]+=1
    q=p**m; print(f"p={p} m={m} (order {2*q}, lambda={2*q//p}, |K|={p**(r-1)}): {len(hist)} weight types")
    for ws,c in sorted(hist.items(),key=lambda t:-t[1]):
        print(f"   {c:6d} hyperplanes  weights {ws}  counting-solvable={count_ok(ws,p**(r-1))}  min weight={min(ws)}")
