"""Check of the M+2m >= |K| lemma on the Jungnickel and quadric families:
every K with max weight + 2*min weight >= |K| must (a) have solvable counting equation and
(b) admit a tiling of K by at most 3 slice translates (a colouring with <= 3 vectors coloured)."""
import itertools, io, contextlib
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from quadric_gh import build
    from quadric_ks import dual_subspaces, count_ok
def tiles_upto3(coords,A,p,r):
    K=[v for v in itertools.product(range(p),repeat=r) if all(sum(a*b for a,b in zip(f,v))%p==0 for f in A)]
    Ks=set(K); groups=defaultdict(list)
    for z in coords: groups[tuple(sum(a*b for a,b in zip(f,z))%p for f in A)].append(z)
    tiles=[frozenset(tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs) for rs in groups.values()]
    add=lambda t,a: tuple((x+y)%p for x,y in zip(t,a))
    full=frozenset(K)
    placed=[frozenset(add(t,a) for a in T) for T in tiles for t in K]
    placed=list(set(placed))
    bysize=sorted(placed,key=len,reverse=True)
    for P in bysize:
        if P==full: return True
    # two or three pieces: first piece must contain 0 (WLOG some piece covers 0)
    zero=tuple([0]*r)
    firsts=[P for P in placed if zero in P]
    for P1 in firsts:
        rest=full-P1
        for P2 in placed:
            if P2<=rest:
                if P2==rest: return True
                rem=rest-P2
                if rem in set(placed): return True
    return False
for p,m,dims in [(3,1,[1,2]),(5,1,[1,2]),(7,1,[1,2]),(3,2,[1,2,3,4])]:
    F,T,coords,_,_=build(p,m); r=len(coords[0])
    for dK in dims:
        n_app=0; bad_count=0; bad_tile=0
        for A in dual_subspaces(p,r,r-dK):
            ws=Counter(tuple(sum(a*b for a,b in zip(f,z))%p for f in A) for z in coords).values()
            M,mn=max(ws),min(ws)
            if M+2*mn>=p**dK:
                n_app+=1
                if not count_ok(list(ws),p**dK): bad_count+=1
                if not tiles_upto3(coords,A,p,r): bad_tile+=1
        print(f"p={p} m={m} dimK={dK}: lemma applies to {n_app} subspaces; counting unsolvable among them: {bad_count}; "
              f"no <=3-piece tiling: {bad_tile}",flush=True)
