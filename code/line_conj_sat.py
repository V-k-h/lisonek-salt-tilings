"""Exhaustive SAT test of the line-slice conjecture for given (p, lambda):
does a normalized GH(Z_p, lambda) exist with a 1-dim direction d in which every slice has 2..p-1 rows?
WLOG d_k = k // lambda (balanced, blocks of size lambda); slice C = {b_C + s d : s in S_C}."""
import sys, itertools, time
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
from collections import Counter

def canon_shape(S,p):  # translation class
    return min(tuple(sorted((x-t)%p for x in S)) for t in range(p))
def shapes_of_size(k,p):
    return sorted({canon_shape(c,p) for c in itertools.combinations(range(p),k)})

def configurations(p,lam):
    """multisets of shapes (sizes 2..p-1) with total size N=p*lam, up to multiplication by units"""
    N=p*lam
    allshapes=[s for k in range(2,p) for s in shapes_of_size(k,p)]
    seen=set(); out=[]
    def rec(i,rem,cur):
        if rem==0:
            key=min(tuple(sorted(canon_shape([(u*x)%p for x in s],p) for s in cur)) for u in range(1,p))
            if key not in seen: seen.add(key); out.append(list(key))
            return
        if i==len(allshapes): return
        s=allshapes[i]
        for c in range(rem//len(s),-1,-1):
            rec(i+1,rem-c*len(s),cur+[s]*c)
    rec(0,N,[]); return out

def equidist_ok(sizes,p,lam):
    """necessary condition (weighted equidistribution): weights split into p groups of total lam (lam>=2)"""
    if lam<2: return True
    if max(sizes)>lam: return False
    sizes=sorted(sizes,reverse=True); bins=[0]*p
    def place(i):
        if i==len(sizes): return all(b==lam for b in bins)
        tried=set()
        for k in range(p):
            if bins[k] in tried or bins[k]+sizes[i]>lam: continue
            tried.add(bins[k]); bins[k]+=sizes[i]
            if place(i+1): return True
            bins[k]-=sizes[i]
        return False
    return place(0)

def build_cnf(p,lam,shapes,zero_slice=0,zero_offset=0):
    N=p*lam; d=[k//lam for k in range(N)]; n=len(shapes)
    vid=[0]
    def new():
        vid[0]+=1; return vid[0]
    X=[[[new() for x in range(p)] for k in range(N)] for C in range(n)]
    cl=[]
    for C in range(n):
        for k in range(N):
            enc=CardEnc.equals(lits=X[C][k],bound=1,top_id=vid[0],encoding=EncType.pairwise); cl+=enc.clauses
    # normalization: the zero row is b_Z + s*d with Z=zero_slice, s=zero_offset, so b_Z = -s*d; column 0 is zero
    for k in range(N): cl.append([X[zero_slice][k][(-zero_offset*d[k])%p]])
    for C in range(n): cl.append([X[C][0][0]])
    for C,D in itertools.combinations(range(n),2):
        E=[[new() for v in range(p)] for k in range(N)]
        for k in range(N):
            for x in range(p):
                for v in range(p):
                    cl.append([-X[C][k][x],-X[D][k][(x-v)%p],E[k][v]])
            enc=CardEnc.equals(lits=E[k],bound=1,top_id=vid[0],encoding=EncType.pairwise); cl+=enc.clauses
        deltas={(s-t)%p for s in shapes[C] for t in shapes[D]}
        for delta in deltas:
            for u in range(p):
                lits=[E[k][(u-delta*d[k])%p] for k in range(N)]
                enc=CardEnc.equals(lits=lits,bound=lam,top_id=vid[0],encoding=EncType.seqcounter)
                vid[0]=max(vid[0],enc.nv); cl+=enc.clauses
        # distinct lines: b_C - b_D is not a multiple of d
        for c in range(p):
            cl.append([-E[k][(c*d[k])%p] for k in range(N)])
    return cl,X,d

def solve(p,lam,shapes,zero_slice=0,zero_offset=0):
    cl,X,d=build_cnf(p,lam,shapes,zero_slice,zero_offset); N=p*lam; n=len(shapes)
    with Cadical153(bootstrap_with=cl) as S:
        ok=S.solve()
        if not ok: return None
        model=set(l for l in S.get_model() if l>0)
    B=[[next(x for x in range(p) if X[C][k][x] in model) for k in range(N)] for C in range(n)]
    rows=[[(B[C][k]+s*d[k])%p for k in range(N)] for C in range(n) for s in shapes[C]]
    return rows

if __name__=="__main__":
    p,lam=int(sys.argv[1]),int(sys.argv[2])
    confs=configurations(p,lam)
    if '--prune' in sys.argv:
        allc=len(confs); confs=[c for c in confs if equidist_ok([len(x) for x in c],p,lam)]
        print(f"pruned by weighted equidistribution: {allc} -> {len(confs)}",flush=True)
    print(f"p={p} lambda={lam} (order {p*lam}): {len(confs)} slice configurations up to symmetry",flush=True)
    found=0; t0=time.time()
    for i,shapes in enumerate(confs):
        t1=time.time(); rows=None
        for z in sorted({shapes.index(sh) for sh in shapes}):
            for off in shapes[z]:
                rows=solve(p,lam,shapes,z,off)
                if rows is not None: break
            if rows is not None: break
        sizes=Counter(len(s) for s in shapes)
        status="UNSAT" if rows is None else "SAT -> candidate counterexample"
        print(f"  [{i+1}/{len(confs)}] sizes {dict(sorted(sizes.items()))} shapes {shapes}: {status} [{time.time()-t1:.1f}s]",flush=True)
        if rows is not None:
            found+=1
            import json; json.dump(rows,open(f"counterexample_p{p}_l{lam}_{i}.json","w"))
    print(f"done: {found} satisfiable configurations [{time.time()-t0:.0f}s]")
