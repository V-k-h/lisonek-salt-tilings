"""lambda=2 line-slice counterexamples in permutation form: which subsets of the orthogonality equations
F_ij(t)=0 (t in {0}, {a_i,-a_j}, {a_i-a_j}) are already unsatisfiable, together with the permutation
constraint (Prop. weighted equidistribution) and the normalizations?
Usage: lam2_subsets.py p config_index families   e.g.  7 3 F0,Fd"""
import sys, itertools, time, threading
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153

def gap_configs(p):
    """multisets of gap classes in F_p^*/{±1} (p elements), up to multiplication by units"""
    classes=sorted({min(a,p-a) for a in range(1,p)})
    rep=lambda a: min(a%p,(-a)%p)
    seen=set(); out=[]
    for combo in itertools.combinations_with_replacement(classes,p):
        key=min(tuple(sorted(rep(u*a) for a in combo)) for u in range(1,p))
        if key not in seen: seen.add(key); out.append(list(key))
    return out

def build(p,gaps,zero_slice,families):
    vid=[0]
    def new(): vid[0]+=1; return vid[0]
    cl=[]
    def exactly(lits,k=1,enc=EncType.seqcounter):
        e=CardEnc.equals(lits=lits,bound=k,top_id=vid[0],encoding=enc); vid[0]=max(vid[0],e.nv); cl.extend(e.clauses)
    Y=[[[[new() for v in range(p)] for eta in range(2)] for b in range(p)] for i in range(p)]   # Y[i][b][eta][v]: entry of row (i,0)
    for i in range(p):
        for b in range(p):
            for eta in range(2): exactly(Y[i][b][eta],enc=EncType.pairwise)
    # normalization: row (zero_slice,0) is zero; column (b=0,eta=0) is zero
    for b in range(p):
        for eta in range(2): cl.append([Y[zero_slice][b][eta][0]])
    for i in range(p): cl.append([Y[i][0][0][0]])
    # permutation constraint: pi_b(i)=Y1-Y0 is a permutation of F_p for every b
    for b in range(p):
        PI=[[new() for v in range(p)] for i in range(p)]
        for i in range(p):
            for x in range(p):
                for v in range(p):
                    cl.append([-Y[i][b][0][x],-Y[i][b][1][(x+v)%p],PI[i][v]])
            exactly(PI[i],enc=EncType.pairwise)
        for v in range(p): exactly([PI[i][v] for i in range(p)])
    for i,j in itertools.combinations(range(p),2):
        E=[[[new() for v in range(p)] for eta in range(2)] for b in range(p)]   # E[b][eta][v]: Y_i - Y_j = v
        for b in range(p):
            for eta in range(2):
                for x in range(p):
                    for v in range(p):
                        cl.append([-Y[i][b][eta][x],-Y[j][b][eta][(x-v)%p],E[b][eta][v]])
                exactly(E[b][eta],enc=EncType.pairwise)
        ts=set()
        if 'F0' in families: ts.add(0)
        if 'Fa' in families: ts|={gaps[i]%p,(-gaps[j])%p}
        if 'Fd' in families: ts.add((gaps[i]-gaps[j])%p)
        star=[int(x[2:]) for x in families if x.startswith('J=')]
        okj = (not star) or (j in star)
        if 'S0' in families and i==0 and okj: ts.add(gaps[0]%p)          # v_0 orthogonal to u_j
        if 'S0b' in families and i==0 and okj: ts.add((-gaps[j])%p)      # u_0 orthogonal to v_j
        for t in ts:
            for u in range(p):
                exactly([E[b][eta][(u-t*b)%p] for b in range(p) for eta in range(2)],k=2)
        # rows of different slices are on different lines: difference of bases not a multiple of d
        for c in range(p):
            cl.append([-E[b][eta][(c*b)%p] for b in range(p) for eta in range(2)])
    return cl

def run(p,gaps,families,limit):
    res=[]
    for z in sorted({gaps.index(g) for g in gaps}):
        S=Cadical153(bootstrap_with=build(p,gaps,z,families))
        timer=threading.Timer(limit,S.interrupt); timer.start()
        r=S.solve_limited(expect_interrupt=True); timer.cancel(); S.delete()
        res.append(r)
        if r is True: return "SAT"
    return "UNSAT" if all(r is False for r in res) else "timeout"

if __name__=="__main__":
    p=int(sys.argv[1]); confs=gap_configs(p)
    idxs=range(len(confs)) if sys.argv[2]=="all" else [int(x) for x in sys.argv[2].split(',')]
    families=sys.argv[3].split(','); limit=float(sys.argv[4]) if len(sys.argv)>4 else 300
    for k in idxs:
        t0=time.time(); r=run(p,confs[k],families,limit)
        print(f"p={p} config {k} gaps {confs[k]} families {families}: {r} [{time.time()-t0:.0f}s]",flush=True)
