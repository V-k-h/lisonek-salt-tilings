import random, sys, itertools
def legendre(a,q): a%=q; return 0 if a==0 else (1 if pow(a,(q-1)//2,q)==1 else -1)
def paley1(q):  # q prime ≡ 3 mod 4, order q+1
    S=[[0]+[1]*q]+[[-1]+[legendre(j-i,q) for j in range(q)] for i in range(q)]
    return [[S[i][j]+(i==j) for j in range(q+1)] for i in range(q+1)]
def paley2(q):  # q prime ≡ 1 mod 4, order 2(q+1)
    Q=[[0]+[1]*q]+[[1]+[legendre(j-i,q) for j in range(q)] for i in range(q)]
    A=[[1,-1],[-1,-1]]; B=[[1,1],[1,-1]]
    n=q+1; H=[[0]*(2*n) for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            blk=B if i==j else ([[Q[i][j]*a for a in r] for r in A])
            for a in range(2):
                for b in range(2): H[2*i+a][2*j+b]=blk[a][b]
    return H
def normalize_bits(H):
    n=len(H); H=[r[:] for r in H]
    for j in range(n):
        if H[0][j]<0:
            for i in range(n): H[i][j]*=-1
    for i in range(n):
        if H[i][0]<0: H[i]=[-v for v in H[i]]
    assert all(sum(H[i][k]*H[j][k] for k in range(n))==(n if i==j else 0) for i in range(n) for j in range(n))
    return [sum(1<<k for k,v in enumerate(r) if v<0) for r in H]
def rank(vs):
    basis=[]
    for v in vs:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    return len(basis)
def dimK(R,parts): return rank([R[a]^R[P[0]] for P in parts for a in P[1:]])
def search(R,q,iters,restarts,seed=0):
    rng=random.Random(seed); n=len(R); best=(10**9,None)
    for _ in range(restarts):
        idx=list(range(n)); rng.shuffle(idx); parts=[idx[i:i+q] for i in range(0,n,q)]
        cur=dimK(R,parts)
        for _ in range(iters):
            g1,g2=rng.sample(range(len(parts)),2); a,b=rng.randrange(q),rng.randrange(q)
            parts[g1][a],parts[g2][b]=parts[g2][b],parts[g1][a]
            d=dimK(R,parts)
            if d<=cur: cur=d
            else: parts[g1][a],parts[g2][b]=parts[g2][b],parts[g1][a]
        if cur<best[0]: best=(cur,[p[:] for p in parts])
    return best
if __name__=="__main__":
  mats={12:("Paley I",paley1(11)),20:("Paley I",paley1(19)),24:("Paley I",paley1(23)),28:("Paley II",paley2(13)),
        36:("Paley II",paley2(17)),44:("Paley I",paley1(43)),48:("Paley I",paley1(47))}
  for n,(name,H) in mats.items():
      R=normalize_bits(H); r=rank(R)
      qs=[q for q in [3,5,7,11] if n%q==0]
      out=[]
      for q in qs:
          d,parts=search(R,q,iters=int(sys.argv[1]),restarts=int(sys.argv[2]))
          generic=n-n//q
          out.append(f"q={q}: best dim K={d} (generic n-n/q={generic}) -> 2^{d} bases, {n//q}x2^{d} rays, saving 2^{r-d}")
      print(f"order {n} ({name}): 2-rank r={r}, Theorem 5 gives 2^{r} bases | "+" | ".join(out) if qs else f"order {n}: no odd prime",flush=True)
