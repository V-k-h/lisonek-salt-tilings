"""Check of the weighted equidistribution identity (1): for K=<d>, columns j,k with d_j=d_k,
sum of slice weights w_C over slices with (r_C)_j-(r_C)_k = t equals lambda for every t."""
import sys
from collections import defaultdict
from p11_ghsearch import normalize
sys.path.insert(0,'quadric')
from quadric_gh import rank_coords
def check(H,p):
    lam=len(H)//p; bad=0; tests=0
    basis,coords=rank_coords(H,p)
    rows=[tuple(r) for r in H]; n=len(rows)
    dirs=set()
    for i in range(n):
        for j in range(i+1,n):
            v=tuple((a-b)%p for a,b in zip(rows[i],rows[j]))
            lead=next(x for x in v if x); inv=pow(lead,-1,p); dirs.add(tuple((x*inv)%p for x in v))
    for d in dirs:
        piv=next(i for i,x in enumerate(d) if x)
        slices=defaultdict(list)
        for r in rows:
            a=(r[piv]*pow(d[piv],-1,p))%p; slices[tuple((x-a*y)%p for x,y in zip(r,d))].append(r)
        if max(len(v) for v in slices.values())<2: continue
        blocks=defaultdict(list)
        for j,x in enumerate(d): blocks[x].append(j)
        for b,J in blocks.items():
            for j in J:
                for k in J:
                    if j>=k: continue
                    tests+=1; tot=defaultdict(int)
                    for base,rs in slices.items(): tot[(rs[0][j]-rs[0][k])%p]+=len(rs)
                    if any(tot[t]!=lam for t in range(p)): bad+=1
    return tests,bad
for k,line in enumerate(open('butson/BH-14-7.txt').read().split()):
    n,p=14,7; dg=line[6:]; H=[[int(dg[i*n+j]) for j in range(n)] for i in range(n)]
    T=B=0
    for i0 in range(n):
        for j0 in range(n):
            t,b=check(normalize(H,p,i0,j0),p); T+=t; B+=b
    print(f"BH(14,7) #{k+1}: {T} (direction, column-pair) tests, {B} violations")
