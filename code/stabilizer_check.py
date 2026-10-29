"""Translation stabilizers of row sets (indexed bases vs distinct contexts)."""
import io, contextlib, sys
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, normalize_bits
    from ks import jung
    sys.path.insert(0,'quadric'); from quadric_gh import build
def stab_bits(R):
    S=set(R); return sum(1 for k in R if k and all((r^k) in S for r in R))
def stab_modp(rows,p):
    S=set(map(tuple,rows)); return sum(1 for k in S if any(k) and all(tuple((a+b)%p for a,b in zip(r,k)) in S for r in S))
for q in [11,23]: print(f"normalized Paley order {q+1}: {stab_bits(normalize_bits(paley1(q)))} nontrivial stabilizer elements")
for p in [3,5,7,11]: print(f"Jungnickel p={p}: {stab_modp(jung(p),p)}")
for p,m in [(3,2),(5,2),(3,3)]:
    F,T,coords,_,_=build(p,m); print(f"quadric q={p}^{m}: {stab_modp(T,p)}")
