import itertools
from collections import Counter
def is_gh(M,p,lam):
    return all(Counter((a-b)%p for a,b in zip(M[i],M[j]))=={e:lam for e in range(p)} for i in range(len(M)) for j in range(i+1,len(M)))
ex5_D=["000000","120201","102210","022011","220110","202101"]
ex5_mod=["000000","012120","021102","022011","001221","010212"]
toM=lambda L:[[int(c) for c in r] for r in L]
print("Example 5 D (as printed) is GH(Z3,2):",is_gh(toM(ex5_D),3,2))
print("Example 5 modified D (as printed) is GH(Z3,2):",is_gh(toM(ex5_mod),3,2))
# formulas with p=3, n=2 (block rows x, cols y in Z_3)
p,n=3,2; i4=pow(4,-1,p); ni=pow(n,-1,p)
def build(c1,c2,c3,c4):
    top=[[ (x*y+c1*x*x)%p for y in range(p)]+[(x*y+c2*x*x)%p for y in range(p)] for x in range(p)]
    bot=[[ (x*y-y*y+c3*x*x)%p for y in range(p)]+[((x*y-y*y+c4*x*x)*ni)%p for y in range(p)] for x in range(p)]
    return top+bot
thm3=build(-i4,n*i4,-i4,-i4)
print("Theorem 3 formulas (p=3,n=2) give:",["".join(map(str,r)) for r in thm3],"GH:",is_gh(thm3,p,2))
cons1=build(0,(n-1)*i4,0,(n-1)*i4)
print("Construction 1 formulas give:",["".join(map(str,r)) for r in cons1],"GH:",is_gh(cons1,p,2))
# which x^2-coefficients (c1..c4) make the bilinear skeleton a GH, for p=5,7?
for p in [5,7,11]:
    sq={(x*x)%p for x in range(1,p)}; n=next(a for a in range(2,p) if a not in sq); ni=pow(n,-1,p); i4=pow(4,-1,p)
    def B(c1,c2,c3,c4):
        return [[ (x*y+c1*x*x)%p for y in range(p)]+[(x*y+c2*x*x)%p for y in range(p)] for x in range(p)]+\
               [[ (x*y-y*y+c3*x*x)%p for y in range(p)]+[((x*y-y*y+c4*x*x)*ni)%p for y in range(p)] for x in range(p)]
    sols=[c for c in itertools.product(range(p),repeat=4) if is_gh(B(*c),p,2)]
    # express solutions as c_i*4 to compare with 'x^2/4' conventions
    print(f"p={p}, n={n}: GH iff 4*(c1,c2,c3,c4) in",sorted({tuple((4*ci)%p for ci in c) for c in sols})[:8],"... total",len(sols),
          "| printed Thm3 is 4c=",((-1)%p,n%p,(-1)%p,(-1)%p))
