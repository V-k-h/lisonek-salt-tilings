import sys, time, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from p7tile import test
    from verify import model, subspaces_dual, wts
p=int(sys.argv[1]); M,_=model(p)
targets=[F for F in subspaces_dual(p,1) if 1 not in wts(M,F,p) and any(w%2 for w in wts(M,F,p))]
print(f"p={p}: {len(targets)} hyperplanes with odd weights and no weight-1 slice",flush=True)
col=0
for i,F in enumerate(targets):
    t=time.time(); ws,ok,use=test(p,F); col+=ok
    print(i,F[0],ws,"colorable" if ok else "NOT COLORABLE (KS!)",dict(use) if use else "",f"{time.time()-t:.1f}s",flush=True)
print(f"DONE p={p}: {col}/{len(targets)} colorable")
