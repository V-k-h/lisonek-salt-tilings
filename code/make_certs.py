"""DRAT certificates for the p=7, lambda=2 line-slice exclusion (Remark on lambda=2).
Usage: make_certs.py INSTANCE_INDEX   (instances = all (configuration, zero-slice, offset) triples)."""
import sys, os, json, time, hashlib, subprocess
from pysat.solvers import Solver
from line_conj_sat import configurations, equidist_ok, build_cnf
p,lam=7,2
DRAT=os.environ.get("DRAT_TRIM")
confs=[c for c in configurations(p,lam) if equidist_ok([len(x) for x in c],p,lam)]
inst=[(ci,z,off) for ci,sh in enumerate(confs) for z in sorted({sh.index(s) for s in sh}) for off in sh[z]]
if sys.argv[1]=="count": print(len(inst)); sys.exit()
ci,z,off=inst[int(sys.argv[1])]; shapes=confs[ci]
tag=f"p7l2_c{ci}_z{z}_o{off}"; cnf=f"certificates/{tag}.cnf"; prf=f"certificates/{tag}.drat"
if os.path.exists(f"certificates/{tag}.json"): print(tag,"already done"); sys.exit()
cl,_,_=build_cnf(p,lam,shapes,z,off)
nv=max(abs(l) for c in cl for l in c)
with open(cnf,"w") as fh:
    fh.write(f"p cnf {nv} {len(cl)}\n")
    for c in cl: fh.write(" ".join(map(str,c))+" 0\n")
t0=time.time(); S=Solver(name="cadical153",with_proof=True,bootstrap_with=cl); sat=S.solve(); ts=time.time()-t0
rec={"tag":tag,"configuration":shapes,"zero_slice":z,"zero_offset":off,"vars":nv,"clauses":len(cl),
     "cnf_sha256":hashlib.sha256(open(cnf,"rb").read()).hexdigest(),"solver":"CaDiCaL 1.5.3 (pysat)","solve_seconds":round(ts,1)}
if sat:
    rec["result"]="SAT"
else:
    with open(prf,"w") as fh:
        for line in S.get_proof(): fh.write(line+"\n")
    rec["proof_sha256"]=hashlib.sha256(open(prf,"rb").read()).hexdigest()
    rec["proof_bytes"]=os.path.getsize(prf)
    t0=time.time(); out=subprocess.run([DRAT,cnf,prf,"-t","200000"],capture_output=True,text=True).stdout
    rec["drat_trim"]=[l for l in out.splitlines() if l.startswith("s ")]; rec["check_seconds"]=round(time.time()-t0,1)
    rec["result"]="UNSAT (verified)" if "s VERIFIED" in out else "UNSAT (proof NOT verified)"
    os.remove(prf)
S.delete()
json.dump(rec,open(f"certificates/{tag}.json","w"),indent=1)
print(tag,rec["result"],rec.get("solve_seconds"),rec.get("check_seconds"),flush=True)
