# Lower-bound case analysis for order 48, checked over every slice-size profile.
# D = dim W_odd + dim H_even <= floor(n_odd/2) + M(E), where E = multiset of even part sizes and
# M(E) = max dim of a binary code on |E| coordinates whose nonzero words have weighted weight >= 12 and = 0 mod 4.
# Need floor(n_odd/2) + M(E) <= 8 for every partition of 48 into parts >= 3.
from functools import lru_cache
N, DMIN = 48, 12

def even_multisets(maxsum, minpart=4):
    yield ()
    for p in range(minpart, maxsum + 1, 2):
        for rest in even_multisets(maxsum - p, p):
            yield (p,) + rest

def M_exact(E, target):
    """True iff some code with all nonzero words allowed has dim > target."""
    n = len(E)
    wt = [sum(E[j] for j in range(n) if v >> j & 1) for v in range(1 << n)]
    ok = [v and wt[v] >= DMIN and wt[v] % 4 == 0 for v in range(1 << n)]
    def dfs(S, last, d):
        if d > target: return True
        for v in range(last + 1, 1 << n):
            if not ok[v]: continue
            if any(v ^ s < v for s in S if s): continue      # v must be the minimum of its coset v+S
            if all(ok[v ^ s] for s in S if s):
                if dfs(S + [v ^ s for s in S], v, d + 1): return True
        return False
    return dfs([0], 0, 0)

worst, checked, needed_exact = 0, 0, []
for E in even_multisets(N):
    s = sum(E)
    if s == N: nodd = 0
    elif N - s >= 6:
        nodd = (N - s) // 3; nodd -= nodd % 2          # odd parts come in even number (N even)
        if len(E) and False: pass
    else: continue
    e = sorted(E)
    # prefix bound: smallest parts with total < DMIN support no codeword
    b, tot = 0, 0
    for x in e:
        if tot + x < DMIN: tot += x; b += 1
        else: break
    pref = len(e) - b
    ham = len(e)
    while ham > 0 and (2 ** ham) * (1 + e.count(4)) > 2 ** len(e): ham -= 1   # sphere bound (radius 5)
    bound = min(pref, ham)
    need = 8 - nodd // 2
    if bound > need:
        needed_exact.append(E)
        if M_exact(tuple(e), need): print("COUNTEREXAMPLE profile", E, "n_odd", nodd); worst = 99
    checked += 1
print("profiles checked:", checked, " needing exact search:", len(needed_exact), needed_exact[:10])
print("all profiles give D <= 8" if worst == 0 else "FAILED")
