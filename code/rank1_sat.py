"""Rank-one lambda=2 counterexamples: H_{(i,e),(b,h)} = x_{ib} + e*a_i*b + h*c_b*i over Z_p
(slices relabelled so pi_b(i) = c_b * i, c_0 = 1). Normalized: row (0,0) and column (0,0) are zero.
Decides whether such a GH(Z_p,2) exists; any solution is verified independently."""
import sys, itertools, time
from collections import Counter
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
p=int(sys.argv[1]); limit=float(sys.argv[2]) if len(sys.argv)>2 else 3600
vid=[0]
def new(): vid[0]+=1; return vid[0]
cl=[]
def exactly(lits,k=1,enc=EncType.seqcounter):
    e=CardEnc.equals(lits=lits,bound=k,top_id=vid[0],encoding=enc); vid[0]=max(vid[0],e.nv); cl.extend(e.clauses)
X=[[[new() for v in range(p)] for b in range(p)] for i in range(p)]
A=[[new() for v in range(p)] for i in range(p)]
C=[[new() for v in range(p)] for b in range(p)]
for i in range(p):
    for b in range(p): exactly(X[i][b],enc=EncType.pairwise)
    exactly(A[i],enc=EncType.pairwise); cl.append([-A[i][0]])
for b in range(p): exactly(C[b],enc=EncType.pairwise); cl.append([-C[b][0]])
cl.append([C[0][1]])
for i in range(p): cl.append([X[i][0][0]])
for b in range(p): cl.append([X[0][b][0]])
# entry variables Y[(i,e)][(b,h)][v]
rows=[(i,e) for i in range(p) for e in range(2)]; cols=[(b,h) for b in range(p) for h in range(2)]
Y={}
for (i,e) in rows:
    for (b,h) in cols:
        Y[(i,e,b,h)]=[new() for v in range(p)]
        exactly(Y[(i,e,b,h)],enc=EncType.pairwise)
        # Y = x_ib + e*a_i*b + h*c_b*i
        for x in range(p):
            for a in range(p):
                for c in range(p):
                    if e==0 and a>0: continue
                    if h==0 and c>0: continue
                    lits=[-X[i][b][x]]
                    if e==1: lits.append(-A[i][a])
                    if h==1: lits.append(-C[b][c])
                    v=(x+(a*b if e else 0)+(c*i if h else 0))%p
                    cl.append(lits+[Y[(i,e,b,h)][v]])
# every pair of rows balanced: each difference value exactly twice over the 2p columns
for r1,r2 in itertools.combinations(rows,2):
    E=[]
    for (b,h) in cols:
        Ev=[new() for v in range(p)]
        for x in range(p):
            for v in range(p):
                cl.append([-Y[r1+(b,h)][x],-Y[r2+(b,h)][(x-v)%p],Ev[v]])
        exactly(Ev,enc=EncType.pairwise); E.append(Ev)
    for v in range(p): exactly([Ev[v] for Ev in E],k=2)
print(f"p={p}: {vid[0]} vars, {len(cl)} clauses",flush=True)
import threading
S=Cadical153(bootstrap_with=cl); tm=threading.Timer(limit,S.interrupt); tm.start(); t0=time.time()
r=S.solve_limited(expect_interrupt=True); tm.cancel()
print(f"p={p}: {'SAT' if r else ('UNSAT' if r is False else 'timeout')} [{time.time()-t0:.0f}s]",flush=True)
if r:
    m=set(l for l in S.get_model() if l>0)
    H=[[next(v for v in range(p) if Y[ro+co][v] in m) for co in cols] for ro in rows]
    ok=all(Counter((a-b)%p for a,b in zip(H[s],H[t]))=={v:2 for v in range(p)} for s in range(2*p) for t in range(s+1,2*p))
    print("verified GH:",ok); print(H)
