"""Sanity check of lemma (*): for every 1-dim K=<d> with a slice of >=2 rows, every slice C with S_C != F_p
has another slice D with S_D - S_C != F_p.  Also check: slice of size p-2 => colourable (dim 1)."""
import itertools, io, contextlib
from collections import defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from quadric_gh import build
for p,m in [(3,1),(5,1),(7,1),(3,2),(5,2)]:
    F,T,coords,_,_=build(p,m); r=len(coords[0]); bad=0; tested=0; sizep2=0; sizep2_ok=0
    seen=set()
    for d in itertools.product(range(p),repeat=r):
        if not any(d): continue
        # normalise d up to scalar
        lead=next(x for x in d if x); inv=pow(lead,-1,p); dn=tuple((x*inv)%p for x in d)
        if dn in seen: continue
        seen.add(dn)
        piv=next(i for i,x in enumerate(dn) if x)
        # slice key: r - a*d with a chosen so coordinate piv becomes 0; S_C records a
        slices=defaultdict(list)
        for z in coords:
            a=z[piv]%p; base=tuple((x-a*y)%p for x,y in zip(z,dn)); slices[base].append(a)
        if max(len(v) for v in slices.values())<2: continue
        tested+=1
        S={k:set(v) for k,v in slices.items()}
        for C,SC in S.items():
            if len(SC)==p: continue
            if not any(len({(b-s)%p for b in SD for s in SC})<p for D,SD in S.items() if D!=C): bad+=1
        sizes=[len(v) for v in S.values()]
        if p-2 in sizes:
            sizep2+=1
            # colourable on Z_p: exact cover of Z_p by translates of -S_C
            tiles=[frozenset((-a)%p for a in v) for v in S.values()]
            def cover(unc):
                if not unc: return True
                x=min(unc)
                for T_ in tiles:
                    for a in T_:
                        P=frozenset((x-a+b)%p for b in T_)
                        if P<=unc and cover(unc-P): return True
                return False
            sizep2_ok+=cover(frozenset(range(p)))
    print(f"p={p} m={m}: 1-dim K with a slice of >=2 rows: {tested}; violations of (*): {bad}; "
          f"with a slice of size p-2: {sizep2}, colourable: {sizep2_ok}")
