import itertools, sys
from collections import defaultdict
from p11_ghsearch import normalize
sys.path.insert(0,'quadric')
from quadric_gh import rank_coords, build
def rank_mod(rows,p):
    rows=[list(r) for r in rows]; rk=0; ncol=len(rows[0]) if rows else 0
    for c in range(ncol):
        piv=next((i for i in range(rk,len(rows)) if rows[i][c]%p),None)
        if piv is None: continue
        rows[rk],rows[piv]=rows[piv],rows[rk]; inv=pow(rows[rk][c],-1,p)
        rows[rk]=[(x*inv)%p for x in rows[rk]]
        for i in range(len(rows)):
            if i!=rk and rows[i][c]%p:
                f=rows[i][c]; rows[i]=[(x-f*y)%p for x,y in zip(rows[i],rows[rk])]
        rk+=1
    return rk
def check(H,p):
    rows=[tuple(r) for r in H]; n=len(rows); N=len(rows[0]); worst=-99; tested=0
    zero=next(i for i,r in enumerate(rows) if not any(r))
    dirs=set()
    for i in range(n):
        for j in range(i+1,n):
            v=tuple((a-b)%p for a,b in zip(rows[i],rows[j])); lead=next(x for x in v if x); inv=pow(lead,-1,p)
            dirs.add(tuple((x*inv)%p for x in v))
    for d in dirs:
        piv=next(i for i,x in enumerate(d) if x); ip=pow(d[piv],-1,p)
        sl=defaultdict(list)
        for idx,r in enumerate(rows):
            a=(r[piv]*ip)%p; sl[tuple((x-a*y)%p for x,y in zip(r,d))].append(idx)
        if max(len(v) for v in sl.values())>=p: continue
        S=list(sl.values()); m=len(S)
        cols=[]
        for j in range(N):
            for k in range(j+1,N):
                if d[j]==d[k]:
                    cols.append([(rows[c[0]][j]-rows[c[0]][k])%p for c in S])
        if not cols: continue
        tested+=1
        r=rank_mod([list(x) for x in zip(*cols)],p)   # rows = slices, columns = column pairs
        worst=max(worst,r-(m//2-1))
    return tested,worst
for k,line in enumerate(open('butson/BH-14-7.txt').read().split()):
    n,p=14,7; dg=line[6:]; H=[[int(dg[i*n+j]) for j in range(n)] for i in range(n)]
    T=0; W=-99
    for i0 in range(0,n,3):
        for j0 in range(0,n,3):
            t,w=check(normalize(H,p,i0,j0),p); T+=t; W=max(W,w)
    print(f"BH(14,7) #{k+1}: {T} directions tested; max of rank - (floor(m/2)-1) = {W} (must be <= 0)")
F,T,coords,_,_=build(5,2)
t,w=check(T,5); print(f"quadric q=25: {t} directions; max excess {w}")
# p=7 permutation enumeration
p=7; u=tuple(range(p))
perms=[(0,)+q for q in itertools.permutations(range(1,p))]
orth=[v for v in perms if sum(a*b for a,b in zip(u,v))%p==0]
mult=[v for v in perms if any(all((c*a-b)%p==0 for a,b in zip(u,v)) for c in range(1,p))]
planes=defaultdict(set)
for v in orth:
    if v in mult: continue
    # canonical plane: set of all F_7-combinations of u and v
    span=frozenset(tuple((x*a+y*b)%p for a,b in zip(u,v)) for x in range(p) for y in range(p))
    planes[span].add(v)
cnt=[sum(1 for w in perms if w in P) for P in planes]
print(f"p=7: perms fixing 0: {len(perms)}, orthogonal to u: {len(orth)}, multiples of u: {len(mult)}, "
      f"others: {len(orth)-len(mult)}, planes: {len(planes)}, permutation vectors per plane: {sorted(set(cnt))}")
