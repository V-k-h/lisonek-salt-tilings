import random, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from modq import paley1, paley2, normalize_bits, rank
def dimK(R,parts): return rank([R[a]^R[P[0]] for P in parts for a in P[1:]])
def slices(R,parts):
    # actual coset partition of rows by K: merge parts whose reps differ by element of K
    vs=[R[a]^R[P[0]] for P in parts for a in P[1:]]; r=rank(vs); out=[]
    for i,P in enumerate(parts):
        for j,Q in enumerate(out):
            if rank(vs+[R[P[0]]^R[parts[Q[0]][0]]])==r: Q.append(i); break
        else: out.append([i])
    return sorted(sum(len(parts[i]) for i in Q) for Q in out)
def search(R,sizes,iters,restarts,seed=0):
    rng=random.Random(seed); n=len(R); best=(10**9,None)
    for _ in range(restarts):
        idx=list(range(n)); rng.shuffle(idx); parts=[]; i=0
        for s in sizes: parts.append(idx[i:i+s]); i+=s
        cur=dimK(R,parts)
        for _ in range(iters):
            g1,g2=rng.sample(range(len(parts)),2); a,b=rng.randrange(len(parts[g1])),rng.randrange(len(parts[g2]))
            parts[g1][a],parts[g2][b]=parts[g2][b],parts[g1][a]
            d=dimK(R,parts)
            if d<=cur: cur=d
            else: parts[g1][a],parts[g2][b]=parts[g2][b],parts[g1][a]
        if cur<best[0]: best=(cur,[p[:] for p in parts])
    return best
for n,H,q,sizelist in [(24,paley1(23),3,[(3,)*8,(6,)*4,(3,3,6,6,6),(12,12),(9,9,6),(6,6,12)]),
                       (48,paley1(47),3,[(3,)*16,(6,)*8,(12,)*4]),
                       (12,paley1(11),3,[(3,)*4])]:
    R=normalize_bits(H)
    for sizes in sizelist:
        d,parts=search(R,sizes,4000,30)
        print(f"order {n}, group sizes {sizes}: best dim K={d}, slice weights {slices(R,parts)}",flush=True)
