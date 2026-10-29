import itertools
from collections import Counter
from ks import jung
from verify import nonsq
def chi(a,p): a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
for p in [3,5,7,11,13]:
    n=nonsq(p); ni=pow(n,-1,p); D=jung(p)
    # paper coords (Thm 7) -> M -> claimed model; check model rows reproduce D via basis b
    b=[D[p-2],D[p-1],D[p+0],D[p+1]]
    i2=pow(2,-1,p); i2n=pow(2*n,-1,p)
    c12=lambda x:((x*x+x)*i2%p, -(x*x+2*x)%p,0,0); c34=lambda y:((y*y-y)*i2n%p,(y-y*y)*ni%p,(1-y)%p,y%p)
    rec=lambda c:[sum(ci*bi[j] for ci,bi in zip(c,b))%p for j in range(2*p)]
    okc=all(rec(c12(x))==D[x] for x in range(p)) and all(rec(c34(y))==D[p+y] for y in range(p))
    M=[[4,1,0,ni],[-2,-1,0,0],[0,0,1,1],[0,0,0,1]]
    Mv=lambda c: tuple(sum(M[i][j]*c[j] for j in range(4))%p for i in range(4))
    okM=all(Mv(c12(x))==(x*x%p,x,0,0) for x in range(p)) and all(Mv(c34(y))==(y*y*ni%p,0,1,y) for y in range(p))
    Z=[(x*x%p,x,0,0) for x in range(p)]+[(y*y*ni%p,0,1,y) for y in range(p)]
    # which subcase do the 2p-1 column hyperplanes fall in? column functional alpha with alpha(Z_i)=D[i][col]
    cases=Counter()
    for col in range(1,2*p) if True else []:
        if col==p: continue  # c^{24}_0 ... keep; skip only c^{13}_0 (col 0)
    cols=[c for c in range(2*p) if c!=0]
    for col in cols:
        sols=[a for a in itertools.product(range(p),repeat=4) if all(sum(ai*zi for ai,zi in zip(a,z))%p==D[i][col] for i,z in enumerate(Z))]
        a=sols[0]
        cases["c13" if col<p else "c24", "alpha1=0" if a[0]==0 else "alpha1!=0"]+=1
    # weight formula check over all alpha with alpha1=1
    okw=True
    for a2,a3,a4 in itertools.product(range(p),repeat=3):
        vals=Counter((z[0]+a2*z[1]+a3*z[2]+a4*z[3])%p for z in Z)
        m1=(-a2*a2*pow(4,-1,p))%p; m2=(a3-n*a4*a4*pow(4,-1,p))%p
        okw&=all(vals.get(v,0)==2+chi(v-m1,p)-chi(v-m2,p) for v in range(p))
    print(p,"paperCoords",okc,"M",okM,"weightFormula",okw,"columnHyperplanes",dict(cases))
