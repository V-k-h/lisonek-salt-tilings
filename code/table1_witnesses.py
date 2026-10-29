"""Witnesses for Table 1: explicit partitions of the rows of normalized Paley matrices into parts of size q,
with dim K_P and slice weights re-verified independently; saved to witnesses/table1_partitions.json."""
import json, io, contextlib
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, paley2, normalize_bits, rank, search
mats=[(12,"Paley I",lambda:paley1(11),3),(20,"Paley I",lambda:paley1(19),5),(24,"Paley I",lambda:paley1(23),3),
      (28,"Paley II",lambda:paley2(13),7),(36,"Paley II",lambda:paley2(17),3),(44,"Paley I",lambda:paley1(43),11),(48,"Paley I",lambda:paley1(47),3)]
import os, sys
ONLY_ORDERS=[int(x) for x in sys.argv[1:]]
out=json.load(open('witnesses/table1_partitions.json')) if os.path.exists('witnesses/table1_partitions.json') else []
done={o['order'] for o in out}
for n,name,mk,q in mats:
    if n in done or (ONLY_ORDERS and n not in ONLY_ORDERS): continue
    R=normalize_bits(mk()); r=rank(R)
    best=None
    for seed in range(4):
        d,parts=search(R,q,iters=4000,restarts=20,seed=seed)
        if best is None or d<best[0]: best=(d,parts)
    d,parts=best
    # independent verification: K = span of within-part differences; slices = cosets of K
    gens=[R[a]^R[P[0]] for P in parts for a in P[1:]]
    basis=[]
    for v in gens:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    K={0}
    for b in basis: K|={x^b for x in K}
    slices=Counter(min(R[i]^k for k in K) for i in range(n))
    ws=sorted(slices.values())
    assert len(basis)==d and all(w%q==0 for w in ws)
    out.append({"order":n,"matrix":name,"q":q,"rank_row_space":r,"dim_K":d,"slice_weights":ws,
                "partition_of_row_indices":[sorted(P) for P in parts],
                "rows_as_bitmasks":R})
    print(f"order {n}: rank {r}, best dim K = {d}, slice weights {ws}",flush=True)
    json.dump(out,open("witnesses/table1_partitions.json","w"),indent=1)
