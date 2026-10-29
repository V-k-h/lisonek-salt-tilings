"""Order-24 minimum: write the 24 normalized Paley rows of order 24 in coordinates of a basis of row(R),
for the exhaustive search h24_min.c over all subspaces K of row(R) of a given codimension.
A KS pair needs every slice of weight >= 3 (Corollary 3.3), and enlarging K to codimension c keeps that,
so if no subspace of codimension 4 (dim 7) has all slice weights >= 3, then 2^8 is the minimum."""
import io, contextlib, sys
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, normalize_bits
R = normalize_bits(paley1(23))
assert len(R) == 24 and R[0] == 0
basis, piv = [], []          # reduced echelon basis of row(R), tracking coordinates
coords_of = {}
def reduce(v):
    c = 0
    for j, (b, p) in enumerate(zip(basis, piv)):
        if (v >> p) & 1:
            v ^= b; c ^= coords_of[j]
    return v, c
for r in R:
    v, c = reduce(r)
    if v:
        p = v.bit_length() - 1
        # keep basis reduced: clear bit p in earlier basis vectors
        j = len(basis)
        coords_of[j] = 0
        basis.append(v); piv.append(p)
        coords_of[j] = 1 << j
        # coordinates of v itself in terms of basis j: we express rows via the original basis vectors below
dim = len(basis)
assert dim == 11, dim
# coordinates of each row: solve r = sum x_j basis_j by elimination on pivots (basis is echelon, pivots distinct)
def coords(r):
    x = 0
    for j, (b, p) in enumerate(zip(basis, piv)):
        if (r >> p) & 1:
            r ^= b; x |= 1 << j
    assert r == 0
    return x
X = [coords(r) for r in R]
# sanity: the map is injective and linear on rows
assert len(set(X)) == 24
for a in range(24):
    for b in range(24):
        ra = R[a] ^ R[b]
        assert coords(ra) == X[a] ^ X[b]
with open(sys.argv[1] if len(sys.argv) > 1 else "h24_rows.txt", "w") as f:
    f.write(" ".join(map(str, X)) + "\n")
print("dim row(R) =", dim, "; wrote", len(X), "row coordinates")
