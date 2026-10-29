import sys, itertools
from collections import Counter
from fq import Fq
FIELDS={(3,1):[0],(5,1):[0],(7,1):[0],(11,1):[0],(11,2):[-1,0],(3,2):[-1,0],(5,2):[2,0],(3,3):[1,1,0],(7,2):[3,0]}  # t^m = ...
def field(p,m):
    mod=[c%p for c in FIELDS[(p,m)]]; F=Fq(p,m,mod)
    # sanity: modpoly irreducible <=> field has no zero divisors
    assert all(F.mul(a,b)!=F.zero for a in F.els if a!=F.zero for b in F.els if b!=F.zero) if p**m<=81 else True
    return F
def jungnickel_Fq(F):
    """LS Construction 1 (normalized Jungnickel GH) over F_q: d1=xy, d2=xy+(n-1)x^2/4, d3=xy-y^2, d4=(xy-y^2+(n-1)x^2/4)/n."""
    n=next(a for a in F.els if not F.is_square(a)); i4=F.inv(F.const(4)); ni=F.inv(n)
    c=F.mul(F.sub(n,F.one),i4); Q=F.els; sq=lambda a: F.mul(a,a); rows=[]
    for x in Q:
        rows.append([F.mul(x,y) for y in Q]+[F.add(F.mul(x,y),F.mul(c,sq(x))) for y in Q])
    for x in Q:
        rows.append([F.sub(F.mul(x,y),sq(y)) for y in Q]+[F.mul(F.add(F.sub(F.mul(x,y),sq(y)),F.mul(c,sq(x))),ni) for y in Q])
    return rows,n
def is_gh(M,diff,elems,lam):
    N=len(M)
    return all(Counter(diff(a,b) for a,b in zip(M[i],M[j]))=={e:lam for e in elems} for i in range(N) for j in range(i+1,N))
def rank_coords(R,p):
    """basis of row space (RREF) and coordinates of each row"""
    basis=[];piv=[]
    for r in R:
        v=list(r)
        for b,c in zip(basis,piv):
            if v[c]: f=v[c]; v=[(x-f*y)%p for x,y in zip(v,b)]
        nz=[i for i,x in enumerate(v) if x]
        if nz:
            c=nz[0]; iv=pow(v[c],-1,p); v=[(x*iv)%p for x in v]
            for k in range(len(basis)):
                if basis[k][c]: f=basis[k][c]; basis[k]=[(x-f*y)%p for x,y in zip(basis[k],v)]
            basis.append(v);piv.append(c)
    coords=[tuple(r[c] for c in piv) for r in R]   # RREF: coordinate = entry at pivot column
    for r,cv in zip(R,coords):
        assert [sum(a*b[j] for a,b in zip(cv,basis))%p for j in range(len(r))]==list(r)
    return basis,coords
def build(p,m):
    F=field(p,m); q=p**m
    MF,n=jungnickel_Fq(F)
    ghF=is_gh(MF,F.sub,F.els,2) if q<=27 else None
    T=[[F.trace(a) for a in row] for row in MF]
    # normalize: subtract first row from all rows, then first column from all columns
    T=[[(a-b)%p for a,b in zip(r,T[0])] for r in T]
    T=[[(a-r[0])%p for a in r] for r in T]
    ghp=is_gh(T,lambda a,b:(a-b)%p,range(p),2*q//p)
    basis,coords=rank_coords(T,p)
    return F,T,coords,ghF,ghp
if __name__=="__main__":
    for p,m in [(3,1),(5,1),(3,2),(5,2),(3,3),(7,2)]:
        F,T,coords,ghF,ghp=build(p,m)
        print(f"p={p} m={m} q={p**m}: GH over F_q: {ghF}; traced matrix GH(Z_{p},{2*p**m//p}): {ghp}; rank={len(coords[0])}",flush=True)
