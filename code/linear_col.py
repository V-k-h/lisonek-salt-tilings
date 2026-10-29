import itertools, sys, io, contextlib
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from verify import model, subspaces_dual, wts
def find(p,alpha,M):
    slices=defaultdict(list)
    for r in M: slices[sum(a*b for a,b in zip(alpha,r))%p].append(r)
    keys=list(slices)
    for beta in itertools.product(range(p),repeat=4):
        B={v:[sum(a*b for a,b in zip(beta,r))%p for r in slices[v]] for v in keys}
        good=[v for v in keys if len(set(B[v]))==len(B[v])]
        for m in range(1,len(good)+1):
            for S in itertools.combinations(good,m):
                if sum(len(B[v]) for v in S)!=p: continue
                # find shifts so translates of B[v] partition Z_p
                sh={}
                def place(i,used):
                    if i==len(S): return True
                    for s in range(p):
                        T={(b+s)%p for b in B[S[i]]}
                        if not T&used:
                            sh[S[i]]=(-s)%p
                            if place(i+1,used|T): return True
                    return False
                if place(0,set()): return beta,list(S),dict(sh),[len(slices[v]) for v in S]
    return None
p=int(sys.argv[1]); M,_=model(p)
targets=[F[0] for F in subspaces_dual(p,1) if 1 not in wts(M,F,p) and any(w%2 for w in wts(M,F,p))]
found=0
from certcheck import check
import time
t0=time.time()
for a in targets:
    res=find(p,a,M); found+=res is not None
    if res and len(sys.argv)>2: assert check(p,a,res[0],set(res[1]),res[2]), ("BAD CERT",a,res)
    if res is None: print("no sliced-linear coloring for",a)
print(f"[{time.time()-t0:.0f}s] p={p}: sliced-linear coloring exists for {found}/{len(targets)}; example:",find(p,targets[0],M))
