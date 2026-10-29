"""q=p^m quadric family: check that every subspace K of row space with dim K <= m+1 is KS-colourable.
Vectorized weight computation over all subspaces (via annihilators); exact-cover test where no slice has weight 1."""
import sys, itertools, time, numpy as np, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from quadric_gh import build
from quadric_ks import colorable
p,m=int(sys.argv[1]),int(sys.argv[2])
ONLY=int(sys.argv[3]) if len(sys.argv)>3 else None
SHARD,NSH=(int(sys.argv[4]),int(sys.argv[5])) if len(sys.argv)>5 else (0,1)
F,T,coords,_,_=build(p,m); r=len(coords[0]); R=np.array(coords,dtype=np.int64)   # rows x r
def annihilators(k):
    """all k x r RREF matrices (row spaces = annihilators of subspaces of dim r-k), yielded in numpy batches"""
    for ci,pivs in enumerate(itertools.combinations(range(r),k)):
        if ci%NSH!=SHARD: continue
        free=[(i,j) for i in range(k) for j in range(r) if j>pivs[i] and j not in pivs]
        nf=len(free)
        base=np.zeros((k,r),dtype=np.int64)
        for i in range(k): base[i,pivs[i]]=1
        total=p**nf; step=200000
        for start in range(0,total,step):
            idx=np.arange(start,min(total,start+step))
            M=np.broadcast_to(base,(len(idx),k,r)).copy()
            for t,(i,j) in enumerate(free):
                M[:,i,j]=(idx//(p**t))%p
            yield M
for dimK in ([ONLY] if ONLY else range(1,m+2)):
    k=r-dimK; t0=time.time(); total=0; no_w1=[]
    for M in annihilators(k):
        lab=(np.einsum('nkr,sr->nsk',M,R)%p)              # n x rows x k
        code=(lab*(p**np.arange(k))).sum(axis=2)             # n x rows
        s=np.sort(code,axis=1)
        eq_prev=np.concatenate([np.zeros((len(s),1),bool),s[:,1:]==s[:,:-1]],axis=1)
        eq_next=np.concatenate([s[:,1:]==s[:,:-1],np.zeros((len(s),1),bool)],axis=1)
        has_single=(~eq_prev & ~eq_next).any(axis=1)
        total+=len(M)
        for idx in np.nonzero(~has_single)[0]: no_w1.append(M[idx].tolist())
    bad=0
    for A in no_w1:
        if not colorable([tuple(z) for z in coords],A,p,r): bad+=1; print("   NOT colourable:",A,flush=True)
    print(f"shard {SHARD}/{NSH} q={p}^{m} dim K={dimK}: {total} subspaces, {len(no_w1)} without a weight-1 slice, "
          f"{len(no_w1)-bad} of those colourable, {bad} not [{time.time()-t0:.0f}s]",flush=True)
