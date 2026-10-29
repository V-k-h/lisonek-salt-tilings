"""Three-slice obstruction test (lambda=2, full permutation model).
Keeps: permutation condition, global F0, distinct lines, and Fa for the three designated slices
against every slice. Designated slices share the gap class [1] (WLOG gap +1 after scaling d and
swapping rows); all other gaps are free variables.
Case A: designated = {1,2,3}, zero row in slice 0 (gap free).  Case B: designated = {0,1,2},
zero row = row (0,0), so a_0 in {+1,-1}.   Usage: lam2_three.py p A|B limit"""
import sys, itertools, time, threading
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
p=int(sys.argv[1]); case=sys.argv[2]; limit=float(sys.argv[3])
vid=[0]
def new(): vid[0]+=1; return vid[0]
cl=[]
def exactly(lits,k=1,enc=EncType.seqcounter,guard=None):
    e=CardEnc.equals(lits=lits,bound=k,top_id=vid[0],encoding=enc); vid[0]=max(vid[0],e.nv)
    cl.extend([c+[-guard] for c in e.clauses] if guard else e.clauses)
Y=[[[[new() for v in range(p)] for eta in range(2)] for b in range(p)] for i in range(p)]
for i in range(p):
    for b in range(p):
        for eta in range(2): exactly(Y[i][b][eta],enc=EncType.pairwise)
G=[[new() for v in range(p)] for i in range(p)]           # gap a_i
for i in range(p): exactly(G[i],enc=EncType.pairwise); cl.append([-G[i][0]])
K=int(sys.argv[4]) if len(sys.argv)>4 else 3     # number of designated slices (0 = control)
desig=(list(range(1,K+1)) if case=='A' else list(range(K)))
DIFF='diff' in sys.argv
for n_,s in enumerate(desig):
    g=2 if (DIFF and n_==1) else 1          # with 'diff', the second designated slice gets class [2]
    if case=='B' and s==0: cl.append([G[0][g],G[0][(-g)%p]])
    else: cl.append([G[s][g]])
for b in range(p):
    for eta in range(2): cl.append([Y[0][b][eta][0]])       # row (0,0) is zero
for i in range(p): cl.append([Y[i][0][0][0]])                # column (0,0) is zero
for b in range(p):                                           # permutation condition
    PI=[[new() for v in range(p)] for i in range(p)]
    for i in range(p):
        for x in range(p):
            for v in range(p): cl.append([-Y[i][b][0][x],-Y[i][b][1][(x+v)%p],PI[i][v]])
        exactly(PI[i],enc=EncType.pairwise)
    for v in range(p): exactly([PI[i][v] for i in range(p)])
for i,j in itertools.combinations(range(p),2):
    E=[[[new() for v in range(p)] for eta in range(2)] for b in range(p)]
    for b in range(p):
        for eta in range(2):
            for x in range(p):
                for v in range(p): cl.append([-Y[i][b][eta][x],-Y[j][b][eta][(x-v)%p],E[b][eta][v]])
            exactly(E[b][eta],enc=EncType.pairwise)
    lits_t=lambda t,u:[E[b][eta][(u-t*b)%p] for b in range(p) for eta in range(2)]
    for u in range(p): exactly(lits_t(0,u),k=2)                       # F0
    for c in range(p): cl.append([-E[b][eta][(c*b)%p] for b in range(p) for eta in range(2)])  # distinct lines
    if i in desig or j in desig:                                      # Fa: t = a_i and t = -a_j, guarded
        for g in range(1,p):
            for u in range(p):
                exactly(lits_t(g,u),k=2,guard=G[i][g])
                exactly(lits_t((-g)%p,u),k=2,guard=G[j][g])
S=Cadical153(bootstrap_with=cl); tm=threading.Timer(limit,S.interrupt); tm.start(); t0=time.time()
r=S.solve_limited(expect_interrupt=True); tm.cancel()
print(f"p={p} case {case} designated={K}{' (classes [1],[2])' if DIFF else ''}: {'SAT' if r else ('UNSAT' if r is False else 'timeout')} [{time.time()-t0:.0f}s]",flush=True)
