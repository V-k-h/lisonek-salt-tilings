import random, sys, itertools
from collections import Counter
sys.setrecursionlimit(10000)
def rand_hadamard(n,rng):
    # randomized backtracking over normalized rows (first row all +1, first column +1)
    cands=[v for v in range(1<<(n-1)) if bin(v).count("1")==n//2]  # bit j -> column j+1 is -1
    rng.shuffle(cands)
    rows=[0]
    def ok(v): return all(bin(v^r).count("1")==n//2 for r in rows)
    def bt():
        if len(rows)==n: return True
        for v in cands:
            if v>rows[-1] and ok(v):
                rows.append(v)
                if bt(): return True
                rows.pop()
        return False
    # random order makes v>rows[-1] restrictive; use sorted-by-random-key ordering instead
    return rows if bt() else None
def rand_hadamard2(n,rng):
    cands=[v for v in range(1<<(n-1)) if bin(v).count("1")==n//2]
    while True:
        rng.shuffle(cands); rows=[0]
        def bt(start):
            if len(rows)==n: return True
            for i in range(start,len(cands)):
                v=cands[i]
                if all(bin(v^r).count("1")==n//2 for r in rows):
                    rows.append(v)
                    if bt(i+1): return True
                    rows.pop()
                    if len(rows)<4: return False   # restart early for diversity
            return False
        if bt(0): 
            return [[1]+[(-1 if (v>>j)&1 else 1) for j in range(n-1)] for v in rows]
def binary_normalizations(H):
    n=len(H); out=set()
    for r in range(n):
        for c in range(n):
            R=tuple(sorted(sum(1<<j for j in range(n) if H[i][j]*H[r][j]*H[i][c]*H[r][c]<0) for i in range(n)))
            out.add(R)
    return out
def span_basis(vs):
    basis=[]
    for v in vs:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    return basis
def tiles(R,V):
    # exact cover of V by translates of R (WLOG one translate is R itself)
    Rs=set(R); Vl=sorted(V); Vset=set(Vl)
    placements={t:frozenset(r^t for r in R) for t in Vl}
    def solve(unc,depth=0):
        if not unc: return True
        x=min(unc)                      # cover smallest uncovered element
        for r in R:
            t=x^r; P=placements[t]
            if P<=unc and solve(unc-P): return True
        return False
    return solve(frozenset(Vl)-placements[0])
if __name__=="__main__":
    n=int(sys.argv[1]); trials=int(sys.argv[2]); rng=random.Random(int(sys.argv[3]) if len(sys.argv)>3 else 0)
    seen={}
    for t in range(trials):
        H=rand_hadamard2(n,rng)
        for R in binary_normalizations(H):
            if R in seen: continue
            B=span_basis(R); V={0}
            for b in B: V|={x^b for x in V}
            seen[R]=(len(B),tiles(R,V))
    stats=Counter(seen.values())
    print(f"order {n}: {len(seen)} distinct normalized binary row sets; (rank, tiles row space) ->",dict(stats))
