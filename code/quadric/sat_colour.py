"""Decide KS-colourability of B(H,K) for a hyperplane K = ker(alpha) by SAT (exact cover of K by slice tiles),
and verify any colouring found against every basis."""
import itertools, sys, random, time
from collections import Counter, defaultdict
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
from quadric_gh import build
from quadric_ks import dual_subspaces, count_ok
def sat_colour(coords,alpha,p):
    r=len(coords[0]); dot=lambda a,z: sum(x*y for x,y in zip(a,z))%p
    K=[k for k in itertools.product(range(p),repeat=r) if dot(alpha,k)==0]
    idx={k:i for i,k in enumerate(K)}
    slices=defaultdict(list)
    for z in coords: slices[dot(alpha,z)].append(z)
    # variable for each vector z = r_C + t (C slice, t in K): z coloured
    var={}; cover=defaultdict(list)
    for v,rs in slices.items():
        A=[tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs]
        for t in K:
            x=len(var)+1; var[(v,t)]=x
            for a in A: cover[idx[tuple((u+w)%p for u,w in zip(t,a))]].append(x)
    top=len(var); clauses=[]
    t0=time.time()
    for k,lits in cover.items():
        enc=CardEnc.equals(lits=lits,bound=1,top_id=top,encoding=EncType.seqcounter)
        top=max(top,enc.nv); clauses+=enc.clauses
    print(f"   encoding: {len(var)} vars, {len(clauses)} clauses [{time.time()-t0:.0f}s]",flush=True)
    with Cadical153(bootstrap_with=clauses) as S:
        ok=S.solve(); model=set(l for l in S.get_model() if l>0) if ok else None
    if not ok: return False,None
    Z=set()
    for (v,t),x in var.items():
        if x in model: Z.add(tuple((a+b)%p for a,b in zip(slices[v][0],t)))
    # independent check: each basis R+k contains exactly one coloured vector
    good=all(sum(tuple((a+b)%p for a,b in zip(z,k)) in Z for z in coords)==1 for k in K)
    return True,good
if __name__=="__main__":
    p,m,per=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
    F,T,coords,_,_=build(p,m); r=len(coords[0]); dot=lambda a,z: sum(x*y for x,y in zip(a,z))%p
    bytype=defaultdict(list)
    for A in dual_subspaces(p,r,1):
        ws=tuple(sorted(Counter(dot(A[0],z) for z in coords).values()))
        if count_ok(ws,p**(r-1)): bytype[ws].append(A[0])
    rng=random.Random(0)
    want=set(tuple(int(x) for x in t.split('-')) for t in sys.argv[4].split(',')) if len(sys.argv)>4 else None
    for ws,L in bytype.items():
        if want and ws not in want: continue
        for a in rng.sample(L,min(per,len(L))):
            t0=time.time(); col,chk=sat_colour(coords,a,p)
            print(f"weights {ws}: alpha={a} colourable={col} verified={chk} [{time.time()-t0:.0f}s]",flush=True)
