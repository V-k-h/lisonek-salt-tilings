# Try to find partitions with all parts >= 3 and D = dim(G cap S_P) > n/6 (i.e. dim K_P < n/3).
import sys, random
exec(open(__file__.replace('falsify.py', 'h48_check.py')).read().split("def sigma")[0].replace("print(", "(lambda *a, **k: None)("))
def D_of(parts):
    diffs = []
    for p in parts:
        for j in p[1:]:
            diffs.append(sum((((g >> p[0]) ^ (g >> j)) & 1) << t for t, g in enumerate(Gb)))
    return k - rank(diffs)
# PSL_2(q) elements as permutations of P^1 (index 0 = infinity)
def act(M, i):
    a, b, c, d = M
    if i == 0: return 0 if c == 0 else 1 + a * pow(c, q - 2, q) % q
    x = i - 1; den = (c * x + d) % q
    return 0 if den == 0 else 1 + (a * x + b) * pow(den, q - 2, q) % q
reps = {}
for a in range(q):
  for b in range(q):
    for c in range(q):
      for d in range(q):
        if (a * d - b * c) % q != 1: continue
        perm = [act((a, b, c, d), i) for i in range(n)]
        cyc, seen = [], set()
        for i in range(n):
            if i in seen: continue
            o = []; j = i
            while j not in seen: seen.add(j); o.append(j); j = perm[j]
            cyc.append(o)
        typ = tuple(sorted(len(o) for o in cyc))
        if typ not in reps and min(typ) >= 3: reps[typ] = cyc
    if len(reps) > 12: break
for typ, cyc in sorted(reps.items()):
    print("cycle type", typ[:4], "... parts", len(cyc), " D =", D_of(cyc), " dim K_P =", n // 2 - D_of(cyc))
# local search maximizing D over partitions with parts >= 3
random.seed(1); best = 0
for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 40):
    lab = list(range(n)); random.shuffle(lab)
    t = random.choice([n // 6, n // 5, n // 4, n // 3])
    part = [i % t for i in range(n)]; random.shuffle(part)
    def parts_of(pt): return [[i for i in range(n) if pt[i] == c] for c in set(pt)]
    cur = D_of(parts_of(part))
    for it in range(3000):
        i, j = random.randrange(n), random.randrange(n)
        new = part[:]
        if random.random() < 0.5: new[i], new[j] = new[j], new[i]
        else: new[i] = new[j]
        ps = parts_of(new)
        if min(len(p) for p in ps) < 3: continue
        d = D_of(ps)
        if d >= cur or random.random() < 0.02: part, cur = new, d
        best = max(best, cur)
print("local search: best D found =", best, " (theorem says <=", n // 6, "for parts >= 3)")
