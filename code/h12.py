import itertools
from ks import rref_basis
from tile import colorable
q=11; QR={(x*x)%q for x in range(1,q)}
chi=lambda a: 0 if a%q==0 else (1 if a%q in QR else -1)
# Paley I: H = I + S with S skew conference, order 12
S=[[0]+[1]*q]+[[-1]+[chi(j-i) for j in range(q)] for i in range(q)]
H=[[S[i][j]+(1 if i==j else 0) for j in range(12)] for i in range(12)]
assert all(sum(H[i][k]*H[j][k] for k in range(12))==(12 if i==j else 0) for i in range(12) for j in range(12))
# normalize: multiply columns so row0 is all +1, then rows so col0 is +1
for j in range(12):
    if H[0][j]<0:
        for i in range(12): H[i][j]*=-1
for i in range(12):
    if H[i][0]<0: H[i]=[-v for v in H[i]]
R=[[0 if v==1 else 1 for v in r] for r in H]
print("F2-rank of rows:",len(rref_basis(R,2)))
def parts(s):
    if not s: yield []; return
    a=s[0]
    for b,c in itertools.combinations(s[1:],2):
        rest=[x for x in s[1:] if x not in (b,c)]
        for P in parts(rest): yield [(a,b,c)]+P
best=None; hist={}
for P in parts(list(range(12))):
    diffs=[[(x+y)%2 for x,y in zip(R[a],R[b])] for (a,b,c) in P]+[[(x+y)%2 for x,y in zip(R[a],R[c])] for (a,b,c) in P]
    B=rref_basis(diffs,2); hist[len(B)]=hist.get(len(B),0)+1
    if best is None or len(B)<len(best[1]): best=(P,B)
print("dim(K) histogram over triple partitions:",hist)
P,B=best; ones=tuple([1]*12)
from ks import span
print("best partition",P,"dim",len(B),"contains all-ones:",ones in span(B,2))
from ks import weights; print("weights",weights(R,B,2))
