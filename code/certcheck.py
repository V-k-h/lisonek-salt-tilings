# independent check of a sliced-linear coloring certificate by brute force over K
import itertools, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from verify import model
def check(p,alpha,beta,S_vals,shift):
    M,_=model(p); dot=lambda a,z: sum(x*y for x,y in zip(a,z))%p
    K=[k for k in itertools.product(range(p),repeat=4) if dot(alpha,k)==0]
    colored=lambda z: dot(alpha,z) in S_vals and dot(beta,z)==shift[dot(alpha,z)]
    return all(sum(colored(tuple((a+b)%p for a,b in zip(r,k))) for r in M)==1 for k in K)
