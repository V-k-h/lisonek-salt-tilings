import itertools, random
from collections import Counter
from ks import rref_basis, span
exec(open('h12.py').read().split("print(\"F2-rank")[0])  # builds R
print("rank",len(rref_basis(R,2)), "sum of rows zero:", all(sum(r[j] for r in R)%2==0 for j in range(12)))
def info(parts):
    diffs=[[(x+y)%2 for x,y in zip(R[P[0]],R[b])] for P in parts for b in P[1:]]
    B=rref_basis(diffs,2) if diffs else []
    Ks=span(B,2) if B else {tuple([0]*12)}
    lab=lambda r: min(tuple((a+b)%2 for a,b in zip(r,k)) for k in Ks)
    return len(B), sorted(Counter(lab(tuple(r)) for r in R).values())
def matchings(s):
    if not s: yield []; return
    for b in s[1:]:
        rest=[x for x in s[1:] if x!=b]
        for M in matchings(rest): yield [(s[0],b)]+M
print("perfect matchings:",Counter(tuple(map(str,info(M))) for M in matchings(list(range(12)))))
print("one triple partition:",info([(0,1,2),(3,4,5),(6,7,8),(9,10,11)]))
print("(6,6) example:",info([tuple(range(6)),tuple(range(6,12))]), " (3,3,6):",info([(0,1,2),(3,4,5),tuple(range(6,12))]), " (5,7):",info([tuple(range(5)),tuple(range(5,12))]))
