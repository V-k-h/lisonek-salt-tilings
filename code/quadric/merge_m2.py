"""Merge the RESULT lines of a sharded dim_m2_check.py run and print the totals in the format of an unsharded run.
Refuses to report a total unless every shard 0..n-1 is present exactly once (a missing shard means a job died).
Usage: merge_m2.py LOG [LOG ...]"""
import sys, json
from collections import Counter

def main(paths):
    res = {}
    for path in paths:
        for line in open(path):
            if line.startswith("RESULT "):
                r = json.loads(line[7:])
                if r["shard"] in res: sys.exit(f"shard {r['shard']} appears twice (second time in {path})")
                res[r["shard"]] = r
    if not res: sys.exit("no RESULT lines found")
    first = next(iter(res.values())); p, d, nsh = first["p"], first["d"], first["nsh"]
    if any((r["p"], r["d"], r["nsh"]) != (p, d, nsh) for r in res.values()): sys.exit("logs come from different runs")
    missing = sorted(set(range(nsh)) - set(res))
    if missing: sys.exit(f"{len(missing)} of {nsh} shards missing: {missing[:20]}{' ...' if len(missing) > 20 else ''}")
    stats = Counter(); all2 = Counter(); fails = []
    for r in res.values():
        stats.update(r["stats"]); all2.update(r["all2"]); fails += r["fails"]
    secs = [r["secs"] for r in res.values()]
    print(f"q={p}^2 dim K={d}: {sum(r['total'] for r in res.values())} subspaces "
          f"[{nsh} shards, longest {max(secs)}s, {sum(secs)} cpu-s in all]")
    for key, n in sorted(stats.items()): print(f"   {key}: {n}")
    print("   all slice weights 2:", {tuple(key.split(" | ")): n for key, n in all2.items()})
    print("   counterexamples to Conjecture 2.2:", len(fails))
    for A, ws in fails: print("   COUNTEREXAMPLE:", A, ws)

if __name__ == "__main__":
    main(sys.argv[1:])
