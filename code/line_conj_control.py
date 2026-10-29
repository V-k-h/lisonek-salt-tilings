"""Positive control for line_conj_sat: allowing singleton slices, the encoding must find the known
BH(10,5) slice profiles 1^8 2 and 1^6 2^2, and the returned matrices must be GH(Z_5,2)."""
from line_conj_sat import solve
from p11_ghsearch import is_gh
p,lam=5,2
for shapes in [[(0,)]*8+[(0,1)], [(0,)]*6+[(0,1),(0,1)], [(0,)]*6+[(0,1),(0,2)], [(0,)]*4+[(0,1)]*3]:
    rows=None
    for z in sorted({shapes.index(s) for s in shapes}):
        for off in shapes[z]:
            rows=solve(p,lam,shapes,z,off)
            if rows: break
        if rows: break
    print(shapes, "SAT" if rows else "UNSAT", "| GH verified:" if rows else "", is_gh(rows,p) if rows else "")
