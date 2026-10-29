"""Conjecture for every 1-dim shift space K=<d> of row(Tr D_q): exact tiling of Z_p by translates of -S_C."""
import sys, itertools, io, contextlib, time
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from quadric_gh import build
def tiles_Zp(shapes,p):
    shapes=[frozenset(s) for s in shapes]
    def cover(unc):
        if not unc: return True
        x=min(unc)
        for S in shapes:
            for a in S:
                P=frozenset((x-a+b)%p for b in S)
                if P<=unc and cover(unc-P): return True
        return False
    return cover(frozenset(range(p)))
def count_ok(ws,N):
    ok=[True]+[False]*N
    for s in range(1,N+1): ok[s]=any(s>=w and ok[s-w] for w in set(ws))
    return ok[N]
p,m=int(sys.argv[1]),int(sys.argv[2])
t0=time.time(); F,T,coords,_,_=build(p,m); r=len(coords[0])
print(f"p={p} m={m}: {len(coords)} rows in F_{p}^{r} [{time.time()-t0:.0f}s]",flush=True)
stats=Counter(); fails=[]; t0=time.time()
for piv in range(r):
    for tail in itertools.product(range(p),repeat=r-piv-1):
        d=(0,)*piv+(1,)+tail
        slices=defaultdict(set)
        for z in coords:
            a=z[piv]; slices[tuple((x-a*y)%p for x,y in zip(z,d))].add(a)
        sizes=tuple(sorted(len(v) for v in slices.values()))
        if 1 in sizes or p in sizes: stats["size 1 or p"]+=1; continue
        if not count_ok(sizes,p): stats["counting unsolvable (KS)"]+=1; continue
        ok=tiles_Zp([{(-a)%p for a in v} for v in slices.values()],p)
        stats[("solvable, nontrivial, colourable="+str(ok),)]+=1
        if not ok: fails.append((d,sizes)); print("  COUNTEREXAMPLE",d,sizes,flush=True)
print(dict(stats),f"[{time.time()-t0:.0f}s]")
