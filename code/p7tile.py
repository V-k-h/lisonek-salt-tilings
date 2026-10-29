import itertools, sys
from collections import Counter
from verify import model, subspaces_dual
def algx(X,Y,sol):
    if not X: yield list(sol); return
    c=min(X,key=lambda c:len(X[c]))
    for r in list(X[c]):
        sol.append(r); cols=select(X,Y,r)
        yield from algx(X,Y,sol)
        deselect(X,Y,r,cols); sol.pop()
def select(X,Y,r):
    cols=[]
    for j in Y[r]:
        for i in X[j]:
            for k in Y[i]:
                if k!=j: X[k].remove(i)
        cols.append(X.pop(j))
    return cols
def deselect(X,Y,r,cols):
    for j in reversed(Y[r]):
        X[j]=cols.pop()
        for i in X[j]:
            for k in Y[i]:
                if k!=j: X[k].add(i)
def test(p,F):
    M,_=model(p); f=F[0]
    piv=next(i for i in range(4) if f[i]); 
    # K = ker f in F_p^4; coordinates: drop pivot coordinate
    K=[v for v in itertools.product(range(p),repeat=4) if sum(a*b for a,b in zip(f,v))%p==0]
    idx={k:i for i,k in enumerate(K)}
    groups={}
    for r in M: groups.setdefault(sum(a*b for a,b in zip(f,r))%p,[]).append(r)
    Y={}
    for g,rs in groups.items():
        A=[tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs]
        for t in K: Y[(g,t)]=[idx[tuple((a+b)%p for a,b in zip(t,s))] for s in A]
    X={j:set() for j in range(len(K))}
    for r,cols in Y.items():
        for j in cols: X[j].add(r)
    sol=next(algx(X,Y,[]),None)
    return sorted(len(v) for v in groups.values()), sol is not None, (Counter(len(groups[g]) for g,_ in sol) if sol else None)
if __name__=="__main__":
  p=int(sys.argv[1]); target=tuple(map(int,sys.argv[2].split(',')))
  from verify import wts
  done=0
  for F in subspaces_dual(p,1):
    M,_=model(p)
    if wts(M,F,p)==target:
        print(F, test(p,F), flush=True); done+=1
        if done>=int(sys.argv[3]): break
