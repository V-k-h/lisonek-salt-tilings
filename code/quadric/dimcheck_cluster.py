#!/usr/bin/env python3
"""Exact-cover check of Conjecture 2.2 on every subspace K of dimension d of the row space of the quadric family,
m = 2 (q = p^2), split into shards for a cluster; used for q = 25, d = 4 in the paper (Corollary 8.13).
Trimmed from the project's cluster script; the code below is unchanged.

  dimcheck p d [--nshards N] [--shard i] [--workers W] [--out DIR] [--method sat|algx]
      The subspaces are split round-robin into N shards and each finished shard is saved to DIR as JSON, so a rerun
      skips finished shards.  With --shard, one shard runs; without it, every unfinished shard runs in a pool of W
      processes.  Exact covers are decided by SAT (default) or by Algorithm X.
  merge p d [--out DIR]
      Combine the shards and print the summary; refuses if a shard is missing, and checks the total against the
      Gaussian binomial.
  selftest
      dimcheck and merge for q = 9, d = 4 with both methods (known: 11011 subspaces, 648 KS, 0 counterexamples).

The paper's run: dimcheck 5 4 --nshards 100, then merge 5 4 (log quadric/logs/m2_dim4_q25.log)."""
import sys, os, json, time, itertools, argparse, multiprocessing, glob, subprocess, tempfile
from collections import Counter

def ncpu(default):
    return int(os.environ.get("SLURM_CPUS_PER_TASK", default))

# ======================= finite fields and the forms Q1, Q2 (fq.py, quadric_gh.py, cross_tiling.py) =======================
class Fq:
    """Finite field F_{p^m} as polynomials mod an irreducible; elements are tuples of length m."""
    def __init__(self,p,m,modpoly):   # modpoly: coefficients c0..c_{m-1} with t^m = sum c_i t^i
        self.p,self.m,self.mod=p,m,modpoly
        self.els=list(itertools.product(range(p),repeat=m))
        self.zero=tuple([0]*m); self.one=tuple([1]+[0]*(m-1))
        self._mul={}
    def add(s,a,b): return tuple((x+y)%s.p for x,y in zip(a,b))
    def neg(s,a): return tuple((-x)%s.p for x in a)
    def sub(s,a,b): return s.add(a,s.neg(b))
    def mul(s,a,b):
        key=(a,b)
        if key in s._mul: return s._mul[key]
        p,m=s.p,s.m; prod=[0]*(2*m-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b): prod[i+j]=(prod[i+j]+x*y)%p
        for k in range(2*m-2,m-1,-1):
            c=prod[k]
            if c:
                prod[k]=0
                for i,ci in enumerate(s.mod): prod[k-m+i]=(prod[k-m+i]+c*ci)%p
        r=tuple(prod[:m]); s._mul[key]=r; return r
    def pow(s,a,e):
        r=s.one
        while e:
            if e&1: r=s.mul(r,a)
            a=s.mul(a,a); e>>=1
        return r
    def inv(s,a): return s.pow(a,s.p**s.m-2)
    def is_square(s,a): return a==s.zero or s.pow(a,(s.p**s.m-1)//2)==s.one
    def trace(s,a):
        t=s.zero; x=a
        for _ in range(s.m): t=s.add(t,x); x=s.pow(x,s.p)
        assert all(v==0 for v in t[1:]); return t[0]
    def const(s,k): return tuple([k%s.p]+[0]*(s.m-1))

FIELDS={(3,1):[0],(5,1):[0],(7,1):[0],(11,1):[0],(11,2):[-1,0],(3,2):[-1,0],(5,2):[2,0],(3,3):[1,1,0],(7,2):[3,0]}  # t^m = ...
def field(p,m):
    mod=[c%p for c in FIELDS[(p,m)]]; F=Fq(p,m,mod)
    # sanity: modpoly irreducible <=> field has no zero divisors
    assert all(F.mul(a,b)!=F.zero for a in F.els if a!=F.zero for b in F.els if b!=F.zero) if p**m<=81 else True
    return F

def forms(p, m):
    F = field(p, m)
    n = next(a for a in F.els if not F.is_square(a))
    c = F.mul(F.sub(n, F.one), F.inv(F.const(4))); ni = F.inv(n)
    sq = lambda x: F.mul(x, x)
    Q1 = {x: F.trace(F.mul(c, sq(x))) for x in F.els}
    Q2 = {y: F.trace(F.mul(F.mul(c, ni), sq(y))) for y in F.els}
    return F, Q1, Q2

# ======================= dimcheck: Conjecture 2.2 on all subspaces of dimension d, m = 2 (dim_m2_check.py) =======================
def count_ok(ws,N):
    ok=bytearray(N+1); ok[0]=1
    for s in range(1,N+1):
        for w in set(ws):
            if s>=w and ok[s-w]: ok[s]=1; break
    return bool(ok[N])
def algx(X,Y,sol):
    if not X: return True
    c=min(X,key=lambda c:len(X[c]))
    for r in list(X[c]):
        cols=[]
        for j in Y[r]:
            for i in X[j]:
                for k in Y[i]:
                    if k!=j: X[k].discard(i)
            cols.append(X.pop(j))
        if algx(X,Y,sol): return True
        for j in reversed(Y[r]):
            X[j]=cols.pop()
            for i in X[j]:
                for k in Y[i]:
                    if k!=j: X[k].add(i)
    return False
def colorable(coords,A,p,r):
    K=[v for v in itertools.product(range(p),repeat=r) if all(sum(a*b for a,b in zip(f,v))%p==0 for f in A)]
    idx={k:i for i,k in enumerate(K)}
    groups={}
    for z in coords: groups.setdefault(tuple(sum(a*b for a,b in zip(f,z))%p for f in A),[]).append(z)
    Y={}
    for g,rs in groups.items():
        T=[tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs]
        for t in K: Y[(g,t)]=[idx[tuple((a+b)%p for a,b in zip(t,s))] for s in T]
    X={j:set() for j in range(len(K))}
    for key,cols in Y.items():
        for j in cols: X[j].add(key)
    return algx(X,Y,[])
def colorable_sat(coords,A,p,r):
    """the exact cover of colorable(), decided by CaDiCaL; a model is checked to be an exact cover before answering True.
    Algorithm X is exponential on some q = 25 subspaces (minutes or more), CaDiCaL takes about 0.1 s."""
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Cadical153
    K=[v for v in itertools.product(range(p),repeat=r) if all(sum(a*b for a,b in zip(f,v))%p==0 for f in A)]
    idx={k:i for i,k in enumerate(K)}
    groups={}
    for z in coords: groups.setdefault(tuple(sum(a*b for a,b in zip(f,z))%p for f in A),[]).append(z)
    pls=[]; cover=[[] for _ in K]
    for g,rs in groups.items():
        T=[tuple((a-b)%p for a,b in zip(rs[0],h)) for h in rs]
        for t in K:
            cells=[idx[tuple((a+b)%p for a,b in zip(t,s))] for s in T]; pls.append(cells)
            for c in cells: cover[c].append(len(pls))
    top=len(pls); cls=[]
    for lits in cover:
        enc=CardEnc.equals(lits=lits,bound=1,top_id=top,encoding=EncType.seqcounter); top=max(top,enc.nv); cls+=enc.clauses
    with Cadical153(bootstrap_with=cls) as S:
        if not S.solve(): return False
        mod=S.get_model()
    cnt=[0]*len(K)
    for i,cells in enumerate(pls):
        if mod[i]>0:
            for c in cells: cnt[c]+=1
    assert all(c==1 for c in cnt), "SAT model is not an exact cover"
    return True

def shard_file(outdir, p, d, nshards, shard):
    return os.path.join(outdir, f"dim_q{p*p}_d{d}_shard{shard:04d}of{nshards}.json")

def dim_shard(p, d, nshards, shard, outdir, method="sat"):
    """dim_m2_check.main restricted to the subspaces with index = shard (mod nshards); result saved as JSON."""
    import numpy as np
    sys.setrecursionlimit(10000)
    m = 2; r = 2 * m + 2; k = r - d
    F, Q1, Q2 = forms(p, m)
    rows = [(Q1[x],) + tuple(x) + (0,) + (0,) * m for x in F.els] + [(Q2[y],) + (0,) * m + (1,) + tuple(y) for y in F.els]
    Rm = np.array(rows, dtype=np.int64)
    stats = Counter(); all2 = Counter(); fails = []; t0 = time.time(); total = 0; index = -1; slow = []
    for pivs in itertools.combinations(range(r), k):
        free = [(i, j) for i in range(k) for j in range(r) if j > pivs[i] and j not in pivs]
        for vals in itertools.product(range(p), repeat=len(free)):
            index += 1
            if index % nshards != shard: continue
            A = np.zeros((k, r), dtype=np.int64)
            for i in range(k): A[i, pivs[i]] = 1
            for (i, j), v in zip(free, vals): A[i, j] = v
            total += 1
            lab = (Rm @ A.T) % p
            ws = sorted(Counter(map(tuple, lab)).values())
            e1in = not A[:, 0].any()
            z3sub = bool(np.linalg.matrix_rank(np.vstack([A, np.eye(r, dtype=np.int64)[m + 1]])) == k)  # K in {z3=0}?
            if set(ws) == {2}:
                all2[("e1 in K" if e1in else "e1 not in K", "K in z3=0" if z3sub else "K not in z3=0")] += 1
            if not count_ok(ws, p ** d):
                stats["KS by counting"] += 1; continue
            if 1 in ws:
                stats["solvable, weight-1 slice"] += 1; continue
            t1 = time.time()
            col = (colorable_sat if method == "sat" else colorable)([tuple(z) for z in rows], A.tolist(), p, r)
            if time.time() - t1 > 60: slow.append((A.tolist(), round(time.time() - t1)))
            stats[f"solvable, no weight 1, colourable={col}"] += 1
            if not col:
                fails.append((A.tolist(), ws)); print("   COUNTEREXAMPLE:", A.tolist(), ws, flush=True)
    res = {"p": p, "d": d, "nshards": nshards, "shard": shard, "method": method, "total": total, "stats": dict(stats),
           "all2": [[list(key), n] for key, n in all2.items()], "fails": fails, "slow": slow, "sec": time.time() - t0}
    fn = shard_file(outdir, p, d, nshards, shard)
    with open(fn + ".tmp", "w") as f: json.dump(res, f)
    os.replace(fn + ".tmp", fn)
    print(f"shard {shard}/{nshards}: {total} subspaces, {len(fails)} counterexamples [{res['sec']:.0f}s]", flush=True)
    return shard

def _dim_shard_star(args): return dim_shard(*args)

def gauss_binom(n, k, p):
    num = den = 1
    for i in range(k): num *= p ** (n - i) - 1; den *= p ** (i + 1) - 1
    return num // den

def cmd_dimcheck(A):
    os.makedirs(A.out, exist_ok=True)
    shards = [A.shard] if A.shard is not None else range(A.nshards)
    todo = [s for s in shards if not os.path.exists(shard_file(A.out, A.p, A.d, A.nshards, s))]
    print(f"q={A.p}^2 dim K={A.d}: {A.nshards} shards, running {len(todo)} ({len(shards) - len(todo)} already done)", flush=True)
    if len(todo) == 1: dim_shard(A.p, A.d, A.nshards, todo[0], A.out, A.method)
    elif todo:
        with multiprocessing.Pool(A.workers) as pool:
            for _ in pool.imap_unordered(_dim_shard_star, [(A.p, A.d, A.nshards, s, A.out, A.method) for s in todo]): pass

def cmd_merge(A):
    files = sorted(glob.glob(os.path.join(A.out, f"dim_q{A.p*A.p}_d{A.d}_shard*of*.json")))
    if not files: sys.exit(f"no shard files in {A.out}")
    R = [json.load(open(f)) for f in files]
    N = {x["nshards"] for x in R}
    if len(N) != 1: sys.exit(f"shard files with different shard counts {sorted(N)}; keep one set")
    N = N.pop(); have = {x["shard"] for x in R}; missing = sorted(set(range(N)) - have)
    if missing: sys.exit(f"{len(missing)} of {N} shards missing, e.g. {missing[:10]}")
    stats = Counter(); all2 = Counter(); fails = []; slow = []; total = 0; sec = 0
    for x in R:
        total += x["total"]; stats.update(x["stats"]); fails += x["fails"]; slow += x["slow"]; sec += x["sec"]
        for key, n in x["all2"]: all2[tuple(key)] += n
    expect = gauss_binom(6, A.d, A.p)
    print(f"q={A.p}^2 dim K={A.d}: {total} subspaces [{sec:.0f}s CPU over {N} shards]")
    print(f"   expected number of {A.d}-dim subspaces of F_{A.p}^6: {expect} -> {'complete' if total == expect else 'MISMATCH'}")
    for key, n in sorted(stats.items()): print(f"   {key}: {n}")
    print("   all slice weights 2:", dict(sorted(all2.items())))
    for Am, ws in fails: print("   COUNTEREXAMPLE:", Am, ws)
    if slow: print(f"   {len(slow)} colourability checks took over a minute (slowest {max(t for _, t in slow)}s)")
    print("   counterexamples to Conjecture 2.2:", len(fails), flush=True)


def cmd_selftest(A):
    me = os.path.abspath(__file__); ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for meth in ("sat", "algx"):
            d = os.path.join(tmp, meth)
            out = ""
            for args in (["dimcheck", 3, 4, "--nshards", 4, "--workers", 4, "--out", d, "--method", meth], ["merge", 3, 4, "--out", d]):
                r = subprocess.run([sys.executable, "-u", me] + [str(x) for x in args], capture_output=True, text=True)
                out += r.stdout + r.stderr
            good = all(s in out for s in ["11011 subspaces", "-> complete", "KS by counting: 648", "colourable=True: 400",
                                          "weight-1 slice: 9963", "counterexamples to Conjecture 2.2: 0"])
            print(f"  {'PASS' if good else 'FAIL'}  dimcheck/merge q=9 d=4, method {meth}", flush=True)
            if not good: ok = False; print(out[-3000:])
    print("selftest", "PASSED" if ok else "FAILED", flush=True)
    sys.exit(0 if ok else 1)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("dimcheck"); s.add_argument("p", type=int); s.add_argument("d", type=int)
    s.add_argument("--nshards", type=int, default=1); s.add_argument("--shard", type=int)
    s.add_argument("--workers", type=int, default=ncpu(1)); s.add_argument("--out", default="out")
    s.add_argument("--method", choices=["sat", "algx"], default="sat")
    s = sub.add_parser("merge"); s.add_argument("p", type=int); s.add_argument("d", type=int); s.add_argument("--out", default="out")
    sub.add_parser("selftest")
    A = ap.parse_args()
    {"dimcheck": cmd_dimcheck, "merge": cmd_merge, "selftest": cmd_selftest}[A.cmd](A)

if __name__ == "__main__":
    main()
