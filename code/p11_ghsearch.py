"""Search generalized Hadamard matrices over Z_p (p=11 by default) for 1-dim shift spaces K=<d> in which no
slice has size 1 or p (the only 1-dim cases not settled by Theorem 10.5 / Corollary 10.8)."""
import sys, random, itertools, time
from collections import defaultdict, Counter
from ks import jung
sys.path.insert(0,'quadric')
from quadric_gh import rank_coords
def is_gh(M,p):
    lam=len(M)//p
    return all(Counter((a-b)%p for a,b in zip(M[i],M[j]))=={e:lam for e in range(p)} for i in range(len(M)) for j in range(i+1,len(M)))
def normalize(H,p,i0=0,j0=0):
    return [[(H[i][j]-H[i0][j]-H[i][j0]+H[i0][j0])%p for j in range(len(H[0]))] for i in range(len(H))]
def ksum(A,B,p): return [[(a+b)%p for a in ra for b in rb] for ra in A for rb in B]
def transpose(A): return [list(c) for c in zip(*A)]
def fourier(p): return [[(x*y)%p for y in range(p)] for x in range(p)]
def analyse(H,p):
    """return Counter of slice-size patterns over directions spanned by row differences;
    'open' = directions with no slice of size 1 or p."""
    _,coords=rank_coords(H,p)
    r=len(coords[0]); n=len(coords)
    def normdir(v):
        lead=next(x for x in v if x); inv=pow(lead,-1,p); return tuple((x*inv)%p for x in v)
    bydir=defaultdict(list)
    for i in range(n):
        for j in range(i+1,n):
            bydir[normdir(tuple((a-b)%p for a,b in zip(coords[i],coords[j])))].append((i,j))
    stats=Counter(); opens=[]
    for d,prs in bydir.items():
        parent=list(range(n))
        def f(x):
            while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
            return x
        for i,j in prs: parent[f(i)]=f(j)
        sizes=Counter(f(i) for i in range(n)).values()
        mn,mx=min(sizes),max(sizes)
        key="has size 1 or p" if (mn==1 or mx==p) else "OPEN"
        stats[key]+=1
        if key=="OPEN": opens.append((d,tuple(sorted(sizes))))
    return r,stats,opens
if __name__=="__main__":
    p=int(sys.argv[1]) if len(sys.argv)>1 else 11
    J=jung(p); F=fourier(p)
    assert is_gh(J,p)
    tests=[("J (Jungnickel)",J),("J^T",transpose(J)),("F (+) J",ksum(F,J,p)),("J (+) F",ksum(J,F,p)),("J (+) J",ksum(J,J,p))]
    rng=random.Random(1)
    for name,H in tests:
        N=len(H); t0=time.time()
        norms=[(0,0)]+[(rng.randrange(N),rng.randrange(N)) for _ in range(int(sys.argv[2]) if len(sys.argv)>2 else 6)]
        if N<=2*p: norms=[(i,j) for i in range(N) for j in range(N)]
        tot=Counter(); allopen=[]; ranks=Counter()
        for i0,j0 in norms:
            Hn=normalize(H,p,i0,j0)
            r,st,op=analyse(Hn,p); tot+=st; allopen+=op; ranks[r]+=1
        print(f"{name}: order {N}, lambda {N//p}, {len(norms)} normalizations, row-space dims {dict(ranks)}; "
              f"directions: {dict(tot)} [{time.time()-t0:.0f}s]",flush=True)
        for d,s in allopen[:5]: print("    OPEN direction",d,"slice sizes",s,flush=True)
