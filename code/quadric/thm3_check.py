import sys
from collections import Counter
from quadric_gh import field, is_gh
def mat(F,n,forms):
    Q=F.els; rows=[]
    for x in Q: rows.append([forms[0](x,y) for y in Q]+[forms[1](x,y) for y in Q])
    for x in Q: rows.append([forms[2](x,y) for y in Q]+[forms[3](x,y) for y in Q])
    return rows
def variants(F,n):
    m,a,s,inv=F.mul,F.add,F.sub,F.inv
    i4=inv(F.const(4)); ni=inv(n); sq=lambda u:m(u,u)
    X2=lambda x: m(sq(x),i4)        # x^2/4
    return {
     "LS Theorem 3 as printed": [lambda x,y:s(m(x,y),X2(x)), lambda x,y:a(m(x,y),m(n,X2(x))),
                                 lambda x,y:s(s(m(x,y),sq(y)),X2(x)), lambda x,y:m(s(s(m(x,y),sq(y)),X2(x)),ni)],
     "LS Construction 1":       [lambda x,y:m(x,y), lambda x,y:a(m(x,y),m(s(n,F.one),X2(x))),
                                 lambda x,y:s(m(x,y),sq(y)), lambda x,y:m(a(s(m(x,y),sq(y)),m(s(n,F.one),X2(x))),ni)],
     "Thm 3 with -n x^2/4 in d2":[lambda x,y:s(m(x,y),X2(x)), lambda x,y:s(m(x,y),m(n,X2(x))),
                                 lambda x,y:s(s(m(x,y),sq(y)),X2(x)), lambda x,y:m(s(s(m(x,y),sq(y)),X2(x)),ni)],
    }
for p,mm in [(3,1),(5,1),(7,1),(11,1),(3,2),(5,2),(3,3)]:
    F=field(p,mm); n=next(u for u in F.els if not F.is_square(u))
    res={name:is_gh(mat(F,n,fs),F.sub,F.els,2) for name,fs in variants(F,n).items()}
    print(f"q={p}^{mm}:",res,flush=True)
