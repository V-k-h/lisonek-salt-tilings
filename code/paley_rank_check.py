"""Proposition 9.3: F2-ranks of raw binary Paley matrices, and the matrix identity used for Paley II."""
import io, contextlib, numpy as np
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, paley2, rank, legendre
binm=lambda H: np.array([[1 if v<0 else 0 for v in r] for r in H],dtype=np.int64)
rk=lambda M: rank([int(''.join(map(str,row[::-1])),2) for row in M])
X=np.array([[0,1],[1,0]]); Mm=np.array([[0,1],[1,1]]); J2=np.ones((2,2),dtype=np.int64)
for q in [5,13,17,29,37,41,53,61]:
    B=binm(paley2(q)); n=q+1
    N=np.array([[1 if (x!=y and legendre(y-x,q)==-1) else 0 for y in range(q)] for x in range(q)])
    c=np.zeros((n,n),dtype=np.int64); c[1:,1:]=N
    pred=(np.kron(np.eye(n,dtype=np.int64),X)+np.kron(np.ones((n,n),dtype=np.int64),Mm)+np.kron(c,J2))%2
    print(f"Paley II q={q}: identity {np.array_equal(B,pred)}, rank {rk(B)} (claim {2*n})")
for q in [3,11,19,43,59,67,83,7,23,31,47,71,79]:
    print(f"Paley I q={q}: rank {rk(binm(paley1(q)))} (claim {q if q%8==3 else (q+1)//2})")
