import itertools, sys
from collections import Counter
from ks import jung, rref_basis
def nonsq(p):
    sq={(x*x)%p for x in range(1,p)}; return next(a for a in range(2,p) if a not in sq)
def model(p):
    n=nonsq(p); i2=pow(2,-1,p); i2n=pow(2*n,-1,p)
    return [((x*x+x)*i2%p, x%p, 0, 0) for x in range(p)]+[((y*y-y)*i2n%p, 0, 1, y%p) for y in range(p)], n
def check_model(p):
    # verify parabola coords reproduce D's rows: basis E_i = combos of paper basis rows
    D=jung(p); n=nonsq(p)
    top=lambda x: D[x%p]; bot=lambda y: D[p+(y%p)]
    b=[top(-2),top(-1),bot(0),bot(1)]
    comb=lambda c: [sum(ci*bi[j] for ci,bi in zip(c,b))%p for j in range(2*p)]
    E=[comb((1,-2,0,0)),comb((0,-1,0,0)),comb((0,0,1,0)),comb((0,0,-1,1))]
    M,_=model(p)
    rows=[[sum(z[i]*E[i][j] for i in range(4))%p for j in range(2*p)] for z in M]
    return rows==[[v%p for v in r] for r in D]
def subspaces_dual(p,k):  # k-dim subspaces of dual F_p^4, as lists of k functionals (RREF)
    for pivs in itertools.combinations(range(4),k):
        free=[(i,j) for i in range(k) for j in range(4) if j>pivs[i] and j not in pivs]
        for vals in itertools.product(range(p),repeat=len(free)):
            C=[[0]*4 for _ in range(k)]
            for i in range(k): C[i][pivs[i]]=1
            for (i,j),v in zip(free,vals): C[i][j]=v
            yield C
def wts(M,F,p):
    return tuple(sorted(Counter(tuple(sum(f[i]*r[i] for i in range(4))%p for f in F) for r in M).values()))
for p in [3,5,7,11,13]:
    M,n=model(p)
    out=[f"p={p} model_ok={check_model(p)}"]
    for dimK in (1,2):  # annihilator dim = 4-dimK
        bad=[F for F in subspaces_dual(p,4-dimK) if 1 not in wts(M,F,p)]
        out.append(f"dimK={dimK}: #without weight-1 slice={len(bad)}")
    hist=Counter(wts(M,F,p) for F in subspaces_dual(p,1))
    even=sum(v for k,v in hist.items() if all(w%2==0 for w in k))
    no1={k:v for k,v in hist.items() if 1 not in k}
    out.append(f"hyperplanes: all-even={even} (p(2p-1)={p*(2*p-1)}); no-weight-1 types={no1}")
    print("; ".join(out))
