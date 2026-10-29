"""Check (1) the direct E-basis proof of LS Theorem 7 and (2) the explicit colouring
of every class-2b(iii) hyperplane (p = 3 mod 4, chi(m1-m2) = -1) of row(D)."""
import itertools, sys, random
from ks import jung
def nonsq(p):
    sq={(x*x)%p for x in range(1,p)}; return next(a for a in range(2,p) if a not in sq)
def chi(a,p): a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def sqrt_mod(a,p): return next(x for x in range(p) if x*x%p==a%p)

def check_theorem7(p):
    n=nonsq(p); D=jung(p); i4=pow(4,-1,p); ni=pow(n,-1,p)
    # D is a GH(Z_p,2): differences of distinct rows contain each element twice
    gh=all(sorted((a-b)%p for a,b in zip(D[i],D[j]))==sorted(list(range(p))*2)
           for i in range(2*p) for j in range(i+1,2*p))
    E1=[0]*p+[(n-1)*i4%p]*p; E2=[y for y in range(p)]*2
    E3=[(-y*y)%p for y in range(p)]+[(-y*y*ni)%p for y in range(p)]
    E4=[y for y in range(p)]+[(y*ni)%p for y in range(p)]
    comb=lambda c: [sum(ci*e[j] for ci,e in zip(c,(E1,E2,E3,E4)))%p for j in range(2*p)]
    rows_ok=all(comb((x*x,x,0,0))==D[x] for x in range(p)) and all(comb((x*x*ni,0,1,x))==D[p+x] for x in range(p))
    return gh and rows_ok

def colouring_for(p,alpha,n):
    """the construction from the pasted proof; returns (m, c1, v0, Cv0, beta)"""
    _,a2,a3,a4=alpha; i2=pow(2,-1,p); i4=pow(4,-1,p)
    x0=(-a2*i2)%p; m1=(-a2*a2*i4)%p; y0=(-n*a4*i2)%p; m2=(a3-n*a4*a4*i4)%p
    w=lambda v: 2+chi(v-m1,p)-chi(v-m2,p)
    if p==3:
        s=sqrt_mod(n*(m1-m2),p); d=pow(s,-1,p); g=(x0-d*y0)%p
        beta=(0,1,g,d); return m1,0,None,[],beta
    v0=next(v for v in range(p) if v not in (m1,m2) and w(v)==2)
    if chi(v0-m1,p)==1:      # both rows of S_v0 in P12: x = x0 +- t ; pair with S_m1
        t=sqrt_mod(v0-m1,p); s=sqrt_mod(n*(m1-m2),p)
        d=(t*pow(s,-1,p))%p; g=(x0+3*t-d*y0)%p; m=m1; step=(2*t)%p
    else:                    # both rows in P34: y = y0 +- t' ; pair with S_m2
        tp=sqrt_mod(n*(v0-m2),p); r=sqrt_mod(m2-m1,p)
        d=(r*pow(tp,-1,p))%p; g=(x0+3*r-d*y0)%p; m=m2; step=(2*r)%p
    beta=(0,1,g,d)
    return m,step,v0,None,beta

def verify(p,alpha,brute):
    n=nonsq(p); ni=pow(n,-1,p)
    M=[(x*x%p,x,0,0) for x in range(p)]+[(x*x*ni%p,0,1,x) for x in range(p)]
    dot=lambda a,z: sum(u*v for u,v in zip(a,z))%p
    m,step,v0,_,beta=colouring_for(p,alpha,n)
    Bm=sorted(dot(beta,r) for r in M if dot(alpha,r)==m)
    if p==3:
        Cm={0}; Cv=set()
    else:
        Bv=[dot(beta,r) for r in M if dot(alpha,r)==v0]
        # place c - B_m = {0,d,2d} and c - B_v0 = {jd,(j+1)d}, j=3,5,...,p-2 (as in the proof)
        a=next(b for b in Bm if (b+step)%p in Bm and (b+2*step)%p in Bm); Cm={(a+2*step)%p}
        b=next(x for x in Bv if (x+step)%p in Bv); Cv={(b+(j+1)*step)%p for j in range(3,p-1,2)}
    def coloured(z):
        v=dot(alpha,z); bz=dot(beta,z)
        return (v==m and bz in Cm) or (v0 is not None and v==v0 and bz in Cv)
    if brute:   # every shift k in K = ker(alpha)
        K=(k for k in itertools.product(range(p),repeat=4) if dot(alpha,k)==0)
        return all(sum(coloured(tuple((a+b)%p for a,b in zip(r,k))) for r in M)==1 for k in K)
    # value check: beta(r+k)=beta(r)+beta(k), alpha(r+k)=alpha(r); beta(K)=F_p since beta_1=0, alpha_1=1
    for t in range(p):
        cnt=0
        for r in M:
            v=dot(alpha,r); bz=(dot(beta,r)+t)%p
            cnt+=(v==m and bz in Cm) or (v0 is not None and v==v0 and bz in Cv)
        if cnt!=1: return False
    return True

def hyperplanes_2b3(p):
    n=nonsq(p); i4=pow(4,-1,p)
    for a2,a3,a4 in itertools.product(range(p),repeat=3):
        m1=(-a2*a2*i4)%p; m2=(a3-n*a4*a4*i4)%p
        if m1!=m2 and chi(m1-m2,p)==-1: yield (1,a2,a3,a4)

if __name__=="__main__":
    for p in [3,7,11,19,23,31,43,47,59,67,71,79,83]:
        t7=check_theorem7(p) if p<=47 else "skipped"
        H=list(hyperplanes_2b3(p)); brute=p<=11
        ok=sum(verify(p,a,brute) for a in H)
        extra=""
        if p in (19,23):
            sample=random.Random(1).sample(H,15); extra=f"; brute-force sample {sum(verify(p,a,True) for a in sample)}/15"
        print(f"p={p}: Thm7 E-basis+GH ok={t7}; class 2b(iii) hyperplanes={len(H)} (expected {p*p*(p-1)//2}); "
              f"explicit colouring valid for {ok}/{len(H)} ({'all shifts' if brute else 'all beta-values'}){extra}",flush=True)
