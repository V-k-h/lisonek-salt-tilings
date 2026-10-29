exec(open(__file__.replace('sextet.py', 'h48_check.py')).read().split("def sigma")[0].replace("print(", "(lambda *a, **k: None)("))
words = []; w = 0
for i in range(1, 1 << k):
    w ^= Gb[(i & -i).bit_length() - 1]
    if w.bit_count() == 8 and (w & 0b1111) == 0b1111: words.append(w)
tetrads = [0b1111] + [o ^ 0b1111 for o in words]
parts = [[j for j in range(n) if t >> j & 1] for t in tetrads]
print("octads through a fixed 4-set:", len(words), " tetrad sizes:", [len(p) for p in parts], " disjoint cover:", sum(map(len, parts)) == n)
diffs = [sum((((g >> p[0]) ^ (g >> j)) & 1) << t for t, g in enumerate(Gb)) for p in parts for j in p[1:]]
D = k - rank(diffs); print("sextet: D =", D, " dim K_P =", n // 2 - D)
