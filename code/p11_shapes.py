"""p=11, dim K=1, slice weights 3 and 4 only: which families of slice shapes (up to translation) admit a
tiling F_11 = (a+T) ⊔ (b+Q1) ⊔ (c+Q2) with T a triple type and Q1,Q2 quadruple types (Q1=Q2 allowed)?"""
import itertools
from collections import defaultdict
p=11
def canon(S):   # translation class representative: lexicographically least translate
    return min(tuple(sorted((x-t)%p for x in S)) for t in range(p))
T3=sorted({canon(c) for c in itertools.combinations(range(p),3)})
T4=sorted({canon(c) for c in itertools.combinations(range(p),4)})
print(len(T3),"triple types,",len(T4),"quadruple types")
idx4={q:i for i,q in enumerate(T4)}
# completions: for triple t (fixed at its representative), complement (8 pts) = disjoint union of two quad translates
pairs=defaultdict(set)
for t in T3:
    comp=frozenset(range(p))-frozenset(t)
    for Q in itertools.combinations(sorted(comp),4):
        rest=comp-frozenset(Q)
        a,b=idx4[canon(Q)],idx4[canon(rest)]
        pairs[t].add((min(a,b),max(a,b)))
for t in T3:
    loops=sum(1 for a,b in pairs[t] if a==b)
    print(f"triple {t}: {len(pairs[t])} completing quad-type pairs ({loops} using one type twice)")
# symmetry: multipliers u in F_11^* act on types
def act(u,S): return canon([(u*x)%p for x in S])
orbit3={}
for t in T3:
    orbit3.setdefault(min(act(u,t) for u in range(1,p)),[]).append(t)
print("triple types up to affine maps:",len(orbit3))
# For each single triple type, maximum number of quad types usable without any completion (max independent set
# in graph with edges = completing pairs, loops = forbidden single types)
def max_indep(edges,n):
    loops={a for a,b in edges if a==b}
    adj=defaultdict(set)
    for a,b in edges:
        if a!=b: adj[a].add(b); adj[b].add(a)
    verts=[v for v in range(n) if v not in loops]
    best=[]
    def bk(R,P,X):
        nonlocal best
        if not P and not X:
            if len(R)>len(best): best=R[:]
            return
        if len(R)+len(P)<=len(best): return
        for v in list(P):
            bk(R+[v],[w for w in P if w not in adj[v] and w!=v],[w for w in X if w not in adj[v]])
            P=[w for w in P if w!=v]; X=X+[v]
    bk([],verts,[])
    return best
print("\nper triple-orbit: max # of quad types a non-tiling family can contain")
for rep,ts in sorted(orbit3.items()):
    b=max_indep(pairs[rep],len(T4))
    print(f"  {rep} (orbit size {len(ts)}): {len(b)} of {len(T4)} quad types")
# a family with ALL triple types present
allpairs=set().union(*pairs.values())
b=max_indep(allpairs,len(T4))
print(f"\nwith all 15 triple types present: largest non-tiling quad family has {len(b)} types: {[T4[i] for i in b]}")
# quad types that can never be in any non-tiling family containing triple t: loops
for rep in sorted(orbit3):
    loops=[T4[a] for a,bb in pairs[rep] if a==bb]
    print(f"  triple {rep}: quad types that tile its complement alone (2 copies): {loops}")
