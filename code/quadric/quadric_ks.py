import itertools, sys, time
from collections import Counter
from quadric_gh import build
def dual_subspaces(p,r,k):
    for pivs in itertools.combinations(range(r),k):
        free=[(i,j) for i in range(k) for j in range(r) if j>pivs[i] and j not in pivs]
        for vals in itertools.product(range(p),repeat=len(free)):
            C=[[0]*r for _ in range(k)]
            for i in range(k): C[i][pivs[i]]=1
            for (i,j),v in zip(free,vals): C[i][j]=v
            yield C
def count_ok(ws,N):
    ok=bytearray(N+1); ok[0]=1
    for s in range(1,N+1):
        for w in set(ws):
            if s>=w and ok[s-w]: ok[s]=1; break
    return bool(ok[N])
def algx(X,Y,sol):
    if not X: return True
    c=min(X,key=lambda c:len(X[c]))
    for r in list(X[c]):
        cols=[]
        for j in Y[r]:
            for i in X[j]:
                for k in Y[i]:
                    if k!=j: X[k].discard(i)
            cols.append(X.pop(j))
        if algx(X,Y,sol): return True
        for j in reversed(Y[r]):
            X[j]=cols.pop()
            for i in X[j]:
                for k in Y[i]:
                    if k!=j: X[k].add(i)
    return False
def colorable(coords,A,p,r):
    K=[v for v in itertools.product(range(p),repeat=r) if all(sum(a*b for a,b in zip(f,v))%p==0 for f in A)]
    idx={k:i for i,k in enumerate(K)}
    groups={}
    for z in coords: groups.setdefault(tuple(sum(a*b for a,b in zip(f,z))%p for f in A),[]).append(z)
    Y={}
    for g,rs in groups.items():
        T=[tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs]
        for t in K: Y[(g,t)]=[idx[tuple((a+b)%p for a,b in zip(t,s))] for s in T]
    X={j:set() for j in range(len(K))}
    for key,cols in Y.items():
        for j in cols: X[j].add(key)
    return algx(X,Y,[])
if __name__=="__main__":
    p,m=int(sys.argv[1]),int(sys.argv[2]); dims=[int(d) for d in sys.argv[3].split(',')]
    F,T,coords,_,_=build(p,m); r=len(coords[0]); print(f"p={p} m={m}: {len(coords)} rows in F_{p}^{r}",flush=True)
    for dK in dims:
        t0=time.time(); stats=Counter(); fails=[]; ks=[]
        for A in dual_subspaces(p,r,r-dK):
            ws=sorted(Counter(tuple(sum(a*b for a,b in zip(f,z))%p for f in A) for z in coords).values())
            solv=count_ok(ws,p**dK)
            if not solv: stats[("KS (counting unsolvable)",tuple(ws))]+=1; continue
            if 1 in ws: stats["solvable, weight-1 slice"]+=1; continue
            col=colorable(coords,A,p,r)
            stats[("solvable, no weight 1, colorable="+str(col),tuple(ws))]+=1
            if not col: fails.append((A,ws)); print("   COUNTEREXAMPLE:",A,ws,flush=True)
        print(f" dim K={dK} (|K|={p**dK}): {dict(stats)}  [{time.time()-t0:.0f}s]",flush=True)
