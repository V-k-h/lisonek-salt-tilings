"""Thm 8.11, type E2 for q = 49: every E2 subspace among 30000 random 4-dim subspaces without e1 (seed 7) is checked for a
KS colouring by SAT (CaDiCaL, conflict budget per instance; each model checked to be an exact cover).
Usage (from code/quadric): dim4_m2_e2_sample.py [conflict budget, default 3000000]"""
import sys, random, time, itertools, threading
from multiprocessing import Pool
from collections import Counter
sys.path.insert(0, '.'); import dim4_m2_cases as D
p = 7; LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 3000000
S = D.Setup(p)

def col(A):
    from pysat.card import CardEnc, EncType
    from pysat.solvers import Cadical153
    K = [v for v in itertools.product(range(p), repeat=6) if all(sum(a * b for a, b in zip(f, v)) % p == 0 for f in A)]
    idx = {k: i for i, k in enumerate(K)}; groups = {}
    for z in S.rows: groups.setdefault(tuple(sum(a * b for a, b in zip(f, z)) % p for f in A), []).append(z)
    pls = []; cover = [[] for _ in K]
    for g, rs in groups.items():
        T = [tuple((a - b) % p for a, b in zip(rs[0], h)) for h in rs]
        for t in K:
            cells = [idx[tuple((a + b) % p for a, b in zip(t, s))] for s in T]; pls.append(cells)
            for c in cells: cover[c].append(len(pls))
    top = len(pls); cls = []
    for lits in cover:
        enc = CardEnc.equals(lits=lits, bound=1, top_id=top, encoding=EncType.seqcounter); top = max(top, enc.nv); cls += enc.clauses
    t0 = time.time()
    with Cadical153(bootstrap_with=cls) as s:
        s.conf_budget(int(LIMIT))
        r = s.solve_limited()
        if r is True:
            mod = s.get_model(); cnt = [0] * len(K)
            for i, cells in enumerate(pls):
                if mod[i] > 0:
                    for c in cells: cnt[c] += 1
            assert all(c == 1 for c in cnt)
    return {True: "colourable", False: "NOT colourable", None: "budget exhausted"}[r], time.time() - t0

def work(item):
    case, A = item; res, sec = col(A); return case, A, res, sec

if __name__ == "__main__":
    rnd = random.Random(p); items = []; n = 0
    while n < 30000:
        A = [[1] + [rnd.randrange(p) for _ in range(5)], [0] + [rnd.randrange(p) for _ in range(5)]]
        if not any(A[1]): continue
        n += 1; case, pred, _ = D.classify(S, A[0], A[1])
        if pred == "E2": items.append((case, A))
    print(f"{len(items)} E2 subspaces in the sample: {dict(Counter(c for c, _ in items))}", flush=True)
    tally = Counter(); t0 = time.time(); slow = []
    with Pool(8) as pool:
        for i, (case, A, res, sec) in enumerate(pool.imap_unordered(work, items), 1):
            tally[(case, res)] += 1
            if res != "colourable" or sec > 60: slow.append((case, A, res, round(sec)))
            if i % 50 == 0 or i == len(items): print(f"[{time.time()-t0:.0f}s] {i}/{len(items)} {dict(tally)}", flush=True)
    for s in slow: print("  slow/unsolved:", s)
