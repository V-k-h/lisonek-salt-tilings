"""Rank-one lambda=2 model with guarded constraints: which interactions force a contradiction?
Row pairs (i,e),(j,e') with i!=j are 'same' if a_i = +-a_j, else 'diff'; column pairs (b,h),(b',h')
with b!=b' are 'same' if c_b = +-c_b', else 'diff'. Usage: rank1_split.py p ROWS COLS limit
with ROWS, COLS in {all, same, diff, none}."""
import sys, itertools, time, threading
from pysat.card import CardEnc, EncType
from pysat.solvers import Cadical153
p=int(sys.argv[1]); RW,CW=sys.argv[2],sys.argv[3]; limit=float(sys.argv[4])
vid=[0]
def new(): vid[0]+=1; return vid[0]
cl=[]
def exactly(lits,k=1,enc=EncType.seqcounter,guard=None):
    e=CardEnc.equals(lits=lits,bound=k,top_id=vid[0],encoding=enc); vid[0]=max(vid[0],e.nv)
    cl.extend([c+[-guard] for c in e.clauses] if guard else e.clauses)
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
rows=[(i,e) for i in range(p) for e in range(2)]; cols=[(b,h) for b in range(p) for h in range(2)]
Y={}
for (i,e) in rows:
    for (b,h) in cols:
        Y[(i,e,b,h)]=[new() for v in range(p)]; exactly(Y[(i,e,b,h)],enc=EncType.pairwise)
        for x in range(p):
            for a in range(p):
                for c in range(p):
                    if (e==0 and a>0) or (h==0 and c>0): continue
                    lits=[-X[i][b][x]]+([-A[i][a]] if e else [])+([-C[b][c]] if h else [])
                    cl.append(lits+[Y[(i,e,b,h)][(x+(a*b if e else 0)+(c*i if h else 0))%p]])
def same_class(V,s,t):
    """literal g <-> (V_s = +-V_t)"""
    g=new(); terms=[]
    for v in range(1,p):
        for w in (v,(-v)%p):
            y=new(); cl.append([-y,V[s][v]]); cl.append([-y,V[t][w]]); cl.append([y,-V[s][v],-V[t][w]]); terms.append(y)
    cl.append([-g]+terms)
    for y in terms: cl.append([g,-y])
    return g
def balanced(vecs1,vecs2,guard):
    E=[]
    for k in range(len(vecs1)):
        Ev=[new() for v in range(p)]
        for x in range(p):
            for v in range(p): cl.append([-vecs1[k][x],-vecs2[k][(x-v)%p],Ev[v]])
        exactly(Ev,enc=EncType.pairwise); E.append(Ev)
    for v in range(p): exactly([Ev[v] for Ev in E],k=2,guard=guard)
SG={}
for i,j in itertools.combinations(range(p),2): SG[('r',i,j)]=same_class(A,i,j)
for b,c in itertools.combinations(range(p),2): SG[('c',b,c)]=same_class(C,b,c)
for r1,r2 in itertools.combinations(rows,2):
    if r1[0]==r2[0] or RW=='none': continue
    g=SG[('r',r1[0],r2[0])]
    guard=None if RW=='all' else (g if RW=='same' else -g)
    balanced([Y[r1+co] for co in cols],[Y[r2+co] for co in cols],guard)
for c1,c2 in itertools.combinations(cols,2):
    if c1[0]==c2[0] or CW=='none': continue
    g=SG[('c',c1[0],c2[0])]
    guard=None if CW=='all' else (g if CW=='same' else -g)
    balanced([Y[ro+c1] for ro in rows],[Y[ro+c2] for ro in rows],guard)
S=Cadical153(bootstrap_with=cl); tm=threading.Timer(limit,S.interrupt); tm.start(); t0=time.time()
r=S.solve_limited(expect_interrupt=True); tm.cancel()
print(f"p={p} rows={RW} cols={CW}: {'SAT' if r else ('UNSAT' if r is False else 'timeout')} [{time.time()-t0:.0f}s]",flush=True)
