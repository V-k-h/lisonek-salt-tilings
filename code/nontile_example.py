"""Example after the non-tiling proposition: A in F_5^2 determines all 6 directions and does not tile F_5^2."""
import itertools
p=5; A=[(0,0),(1,0),(0,1),(1,1),(2,3)]
dirs=set()
for a,b in itertools.combinations(A,2):
    d=((b[0]-a[0])%p,(b[1]-a[1])%p); inv=pow(d[0],-1,p) if d[0] else pow(d[1],-1,p); dirs.add(((d[0]*inv)%p,(d[1]*inv)%p))
G=[(x,y) for x in range(p) for y in range(p)]
tiles=any(len({((a[0]+t[0])%p,(a[1]+t[1])%p) for a in A for t in T})==p*p for T in itertools.combinations(G,p))
print("directions:",len(dirs),"of",p+1,"| tiles F_5^2:",tiles)
