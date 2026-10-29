import itertools, numpy as np, sys
def jung(p):
    sq={(x*x)%p for x in range(1,p)}; n=next(a for a in range(2,p) if a not in sq)
    i4=pow(4,-1,p); ni=pow(n,-1,p)
    rows=[]
    for x in range(p):
        rows.append([ (x*y)%p for y in range(p)]+[(x*y+(n-1)*x*x*i4)%p for y in range(p)])
    for x in range(p):
        rows.append([ (x*y-y*y)%p for y in range(p)]+[((x*y-y*y+(n-1)*x*x*i4)*ni)%p for y in range(p)])
    return [[int(v) for v in r] for r in rows]
def rref_basis(M,p):
    M=[list(r) for r in M]; B=[]; piv=[]
    for r in M:
        r=[v%p for v in r]
        for b,c in zip(B,piv):
            if r[c]: f=r[c]; r=[(a-f*bb)%p for a,bb in zip(r,b)]
        nz=[i for i,v in enumerate(r) if v]
        if nz:
            c=nz[0]; inv=pow(r[c],-1,p); r=[(v*inv)%p for v in r]
            for j in range(len(B)):
                if B[j][c]: f=B[j][c]; B[j]=[(a-f*bb)%p for a,bb in zip(B[j],r)]
            B.append(r); piv.append(c)
    return B
def span(B,p):
    out=set()
    for co in itertools.product(range(p),repeat=len(B)):
        out.add(tuple(int(sum(c*b[i] for c,b in zip(co,B))%p) for i in range(len(B[0]))))
    return out
def subspaces(B,t,p):
    # all t-dim subspaces of span(B), via RREF coefficient matrices
    k=len(B); seen=set()
    for pivs in itertools.combinations(range(k),t):
        free=[(i,j) for i in range(t) for j in range(k) if j>pivs[i] and j not in pivs]
        for vals in itertools.product(range(p),repeat=len(free)):
            C=[[0]*k for _ in range(t)]
            for i in range(t): C[i][pivs[i]]=1
            for (i,j),v in zip(free,vals): C[i][j]=v
            yield [[sum(C[i][j]*B[j][c] for j in range(k))%p for c in range(len(B[0]))] for i in range(t)]
def weights(R,K,p):
    Ks=span(K,p) if K else {tuple([0]*len(R[0]))}
    cos={}
    for r in map(tuple,R):
        key=min(tuple((a-b)%p for a,b in zip(r,k)) for k in Ks)
        cos[key]=cos.get(key,0)+1
    return sorted(cos.values())
if __name__=="__main__":
    p=int(sys.argv[1]); D=jung(p); B=rref_basis(D,p); print("p",p,"rank",len(B))
    for t in range(1,len(B)):
        good=[K for K in subspaces(B,t,p) if all(w%2==0 for w in weights(D,K,p))]
        print("dim",t,"#K with all weights even:",len(good), (weights(D,good[0],p) if good else ""))
        if good: break
