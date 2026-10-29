import sys
from collections import Counter
from p11_ghsearch import is_gh, normalize, analyse
for path,n,p in [(sys.argv[1],10,5),(sys.argv[2],14,7)]:
    for k,line in enumerate(open(path).read().split()):
        digits=line[6:]; assert len(digits)==n*n and set(digits)<=set("0123456789")
        H=[[int(digits[i*n+j]) for j in range(n)] for i in range(n)]
        assert is_gh(H,p), "not GH"
        tot=Counter(); prof=Counter()
        for i0 in range(n):
            for j0 in range(n):
                r,st,op=analyse(normalize(H,p,i0,j0),p); tot+=st
                assert not op
        print(f"BH({n},{p}) #{k+1}: GH verified; all {n*n} normalizations: {dict(tot)}")
