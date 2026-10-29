"""Finite field F_{p^m} as polynomials mod an irreducible; elements are tuples of length m."""
import itertools
class Fq:
    def __init__(self,p,m,modpoly):   # modpoly: coefficients c0..c_{m-1} with t^m = sum c_i t^i
        self.p,self.m,self.mod=p,m,modpoly
        self.els=list(itertools.product(range(p),repeat=m))
        self.zero=tuple([0]*m); self.one=tuple([1]+[0]*(m-1))
        self._mul={}
    def add(s,a,b): return tuple((x+y)%s.p for x,y in zip(a,b))
    def neg(s,a): return tuple((-x)%s.p for x in a)
    def sub(s,a,b): return s.add(a,s.neg(b))
    def scal(s,k,a): return tuple((k*x)%s.p for x in a)
    def mul(s,a,b):
        key=(a,b)
        if key in s._mul: return s._mul[key]
        p,m=s.p,s.m; prod=[0]*(2*m-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b): prod[i+j]=(prod[i+j]+x*y)%p
        for k in range(2*m-2,m-1,-1):
            c=prod[k]
            if c:
                prod[k]=0
                for i,ci in enumerate(s.mod): prod[k-m+i]=(prod[k-m+i]+c*ci)%p
        r=tuple(prod[:m]); s._mul[key]=r; return r
    def pow(s,a,e):
        r=s.one
        while e:
            if e&1: r=s.mul(r,a)
            a=s.mul(a,a); e>>=1
        return r
    def inv(s,a): return s.pow(a,s.p**s.m-2)
    def is_square(s,a): return a==s.zero or s.pow(a,(s.p**s.m-1)//2)==s.one
    def trace(s,a):
        t=s.zero; x=a
        for _ in range(s.m): t=s.add(t,x); x=s.pow(x,s.p)
        assert all(v==0 for v in t[1:]); return t[0]
    def const(s,k): return tuple([k%s.p]+[0]*(s.m-1))
