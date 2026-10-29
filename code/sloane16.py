import sys, io, contextlib
from collections import Counter
from had16 import binary_normalizations, span_basis, tiles
def load(path):
    import re; rows=[l.strip() for l in open(path) if re.fullmatch(r"[+-]+",l.strip())]
    assert all(set(r)<= set("+-") for r in rows)
    H=[[1 if ch=="+" else -1 for ch in r] for r in rows]; n=len(H)
    assert all(len(r)==n for r in H)
    assert all(sum(H[i][k]*H[j][k] for k in range(n))==(n if i==j else 0) for i in range(n) for j in range(n))
    return H
for f in sys.argv[1:]:
    H=load(f); res=Counter()
    for R in binary_normalizations(H):
        B=span_basis(R); V={0}
        for b in B: V|={x^b for x in V}
        res[(len(B),tiles(R,V))]+=1
    print(f.split('/')[-1],"(rank, tiles) over all normalizations:",dict(res),flush=True)
