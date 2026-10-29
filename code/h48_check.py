# Check the order-48 claims: G = colcode(R) + <1> is a doubly-even self-dual [48,24,12] code,
# and the 16 orbits of x -> -1/(x+1) give dim K_P = 16 with the orbits as slices.
import sys
q = int(sys.argv[1]) if len(sys.argv) > 1 else 47
n = q + 1
sq = {(x * x) % q for x in range(1, q)}
chi = lambda a: 0 if a % q == 0 else (1 if a % q in sq else -1)
# rows/cols: index 0 = infinity, index 1+x = x in F_q
H = [[0] * n for _ in range(n)]
H[0] = [1] * n
for x in range(q):
    H[1 + x][0] = -1
    for y in range(q):
        H[1 + x][1 + y] = (1 if x == y else 0) + chi(y - x)
assert all(sum(H[i][k] * H[j][k] for k in range(n)) == (n if i == j else 0) for i in range(n) for j in range(n))
for i in range(n):            # normalize first column (first row is already all +1)
    if H[i][0] == -1:
        H[i] = [-v for v in H[i]]
rows = [sum(1 << k for k in range(n) if H[i][k] == -1) for i in range(n)]   # R_i as bitmasks over columns
cols = [sum(1 << i for i in range(n) if H[i][k] == -1) for k in range(n)]   # columns as bitmasks over rows

def basis(vecs):
    b = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in b: v ^= b[h]
            else: b[h] = v; break
    return b
def rank(vecs): return len(basis(vecs))
def reduce(b, v):
    while v:
        h = v.bit_length() - 1
        if h not in b: return v
        v ^= b[h]
    return 0

ONE = (1 << n) - 1
print("dim row(R) =", rank(rows))
Gb = list(basis(cols + [ONE]).values())
print("dim G =", len(Gb))
pc = lambda v: v.bit_count()
print("self-orthogonal:", all(pc(a & b) % 2 == 0 for a in Gb for b in Gb), " generators doubly-even:", all(pc(a) % 4 == 0 for a in Gb))
# minimum distance by Gray-code enumeration
k = len(Gb); w = 0; mind = n; cnt = {}
for i in range(1, 1 << k):
    w ^= Gb[(i & -i).bit_length() - 1]
    c = pc(w); cnt[c] = cnt.get(c, 0) + 1
print("weight distribution:", dict(sorted(cnt.items())))

def sigma(x):  # x in {None} u F_q
    if x is None: return 0
    if (x + 1) % q == 0: return None
    return (-pow(x + 1, q - 2, q)) % q
idx = lambda x: 0 if x is None else 1 + x
seen = set(); parts = []
for x in [None] + list(range(q)):
    if idx(x) in seen: continue
    orb = []; y = x
    while idx(y) not in seen:
        seen.add(idx(y)); orb.append(idx(y)); y = sigma(y)
    parts.append(orb)
print("orbit sizes:", sorted(len(p) for p in parts))

def D_of(parts):   # dim(G cap S_P) = 24 - rank of column differences of G within parts
    diffs = []
    for p in parts:
        for j in p[1:]:
            diffs.append(sum(((g >> p[0]) & 1 ^ (g >> j) & 1) << t for t, g in enumerate(Gb)))
    return k - rank(diffs)
KP = basis([rows[p[0]] ^ rows[j] for p in parts for j in p[1:]])
print("dim K_P =", len(KP), "  n/2 - D =", n // 2 - D_of(parts))
# slice partition of K_P
lab = {}
for i in range(n):
    for j in lab:
        if reduce(KP, rows[i] ^ rows[j]) == 0: lab[i] = lab[j]; break
    else: lab[i] = i
sl = {}
for i, l in lab.items(): sl.setdefault(l, []).append(i)
print("slice weights:", sorted(len(s) for s in sl.values()), " slices == orbits:", sorted(map(sorted, sl.values())) == sorted(map(sorted, parts)))
