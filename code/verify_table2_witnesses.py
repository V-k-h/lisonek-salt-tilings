"""Check the stored witnesses for Table 2 (witnesses/table1_partitions.json) without rerunning the search:
rebuild each normalized Paley matrix, compare its rows with the stored ones, and recompute dim K_P and the slice
weights of the stored partition by Gaussian elimination over F_2 (coset representatives by full reduction against
an echelon basis of K_P, so K_P itself is never enumerated)."""
import json, io, contextlib
from collections import Counter
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, paley2, normalize_bits, rank

MATRICES = {12: lambda: paley1(11), 20: lambda: paley1(19), 24: lambda: paley1(23), 28: lambda: paley2(13),
            36: lambda: paley2(17), 44: lambda: paley1(43), 48: lambda: paley1(47)}

def echelon(vectors):
    basis = []
    for v in sorted(vectors, reverse=True):
        for b in basis: v = min(v, v ^ b)
        if v: basis.append(v)
    return basis

def representative(v, basis):
    changed = True
    while changed:
        changed = False
        for b in basis:
            if v ^ b < v: v ^= b; changed = True
    return v

ok_all = True
for w in json.load(open("witnesses/table1_partitions.json")):
    n, q, parts = w["order"], w["q"], w["partition_of_row_indices"]
    R = normalize_bits(MATRICES[n]())
    basis = echelon([R[a] ^ R[P[0]] for P in parts for a in P[1:]])
    weights = sorted(Counter(representative(R[i], basis) for i in range(n)).values())
    ok = (R == w["rows_as_bitmasks"] and len(basis) == w["dim_K"] and weights == w["slice_weights"]
          and sorted(len(P) for P in parts) == [q] * (n // q) and all(x % q == 0 for x in weights))
    ok_all &= ok
    print(f"order {n}: rank {rank(R)}, dim K_P = {len(basis)}, slice weights {dict(Counter(weights))}: {'ok' if ok else 'MISMATCH'}")
print("all witnesses verified" if ok_all else "MISMATCH")
