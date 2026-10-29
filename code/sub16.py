import sys, itertools
from collections import Counter
from had16 import binary_normalizations, span_basis
from sloane16 import load
def all_subspaces(B):
    # enumerate subspaces of span(B) via RREF coefficient matrices over F2
    k=len(B)
    for t in range(1,k+1):
        for pivs in itertools.combinations(range(k),t):
            free=[(i,j) for i in range(t) for j in range(k) if j>pivs[i] and j not in pivs]
            for vals in itertools.product((0,1),repeat=len(free)):
                C=[[0]*k for _ in range(t)]
                for i in range(t): C[i][pivs[i]]=1
                for (i,j),v in zip(free,vals): C[i][j]=v
                yield [ (lambda row: __import__('functools').reduce(lambda a,b:a^b,[B[j] for j in range(k) if row[j]],0))(C[i]) for i in range(t)]
def elems(gens):
    S={0}
    for g in gens: S|={x^g for x in S}
    return S
def count_ok(ws,N):
    ok=[True]+[False]*N
    for s in range(1,N+1): ok[s]=any(s>=w and ok[s-w] for w in set(ws))
    return ok[N]
def colorable(R,Kset):
    Kl=sorted(Kset); groups={}
    for r in R:
        key=min(r^k for k in Kl); groups.setdefault(key,[]).append(r)
    rows={}
    for g,rs in groups.items():
        A=[rs[0]^h for h in rs]
        for t in Kl: rows[(g,t)]=frozenset(t^a for a in A)
    X={x:set() for x in Kl}
    for r,P in rows.items():
        for x in P: X[x].add(r)
    def solve(unc):
        if not unc: return True
        x=min(unc,key=lambda c:len([r for r in X[c] if rows[r]<=unc]))
        for r in X[x]:
            if rows[r]<=unc and solve(unc-rows[r]): return True
        return False
    return solve(frozenset(Kl)), sorted(len(v) for v in groups.values())
if __name__=="__main__":
    for f in sys.argv[1:]:
        H=load(f); Rs=sorted(binary_normalizations(H))[:int(4)]
        stats=Counter()
        for R in Rs:
            B=span_basis(R); seenP=set()
            for gens in all_subspaces(B):
                K=elems(gens)
                part=frozenset(frozenset(r for r in R if (r^r0) in K) for r0 in R)
                if part in seenP: continue
                seenP.add(part)
                ws=sorted(len(p) for p in part)
                if min(ws)<=2 or 4 in ws or not count_ok(ws,len(K)): continue
                # test on the minimal K for this partition
                Kmin=elems(span_basis([a^b for p in part for a in p for b in p]))
                col,_=colorable(R,Kmin)
                stats[(tuple(ws),col)]+=1
                if not col: print("  NON-COLORABLE with solvable counting!",f,ws,"dimK",len(span_basis(list(Kmin))),flush=True)
        print(f.split('/')[-1],"interesting partitions (weights, colorable):",dict(stats),flush=True)
