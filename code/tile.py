import sys, itertools
from ks import *
def colorable(R,K,p):
    Ks=sorted(span(K,p)); idx={k:i for i,k in enumerate(Ks)}; N=len(R[0])
    sub=lambda a,b: tuple((x-y)%p for x,y in zip(a,b)); add=lambda a,b: tuple((x+y)%p for x,y in zip(a,b))
    cos={}
    for r in map(tuple,R):
        key=min(sub(r,k) for k in Ks); cos.setdefault(key,[]).append(r)
    placements=[]
    for rs in cos.values():
        r0=rs[0]; A=[sub(r0,h) for h in rs]
        for t in Ks: placements.append(frozenset(idx[add(t,a)] for a in A))
    placements=list(set(placements))
    cover={i:[P for P in placements if i in P] for i in range(len(Ks))}
    def solve(unc):
        if not unc: return True
        i=min(unc,key=lambda j:sum(1 for P in cover[j] if P<=unc))
        for P in cover[i]:
            if P<=unc and solve(unc-P): return True
        return False
    return solve(frozenset(range(len(Ks)))), sorted(len(v) for v in cos.values())
def counting_ok(ws,size):
    ok=[True]+[False]*size
    for s in range(1,size+1): ok[s]=any(s>=w and ok[s-w] for w in set(ws))
    return ok[size]
if __name__=="__main__":
    p=int(sys.argv[1]); t=int(sys.argv[2]); D=jung(p); B=rref_basis(D,p)
    stats={}
    for K in subspaces(B,t,p):
        col,ws=colorable(D,K,p); c=counting_ok(ws,p**t)
        key=(tuple(ws),c,col); stats[key]=stats.get(key,0)+1
    for k,v in sorted(stats.items()): print("weights",k[0],"count-solvable",k[1],"colorable",k[2],"#",v)
