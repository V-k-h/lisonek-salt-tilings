exec(open('h12b.py').read().split('def matchings')[0].split('print("rank"')[0])
exec('def info'+open('h12b.py').read().split('def info')[1].split('def matchings')[0])
import random
random.seed(1)
def randpart(sizes):
    idx=list(range(12)); random.shuffle(idx); out=[]; i=0
    for s in sizes: out.append(tuple(idx[i:i+s])); i+=s
    return out
for sizes in [(3,3,3,3),(3,3,6),(3,5,4),(5,7),(3,3,2,4),(3,3,3,1,2),(2,2,4,4),(6,6),(3,9),(1,11)]:
    res=Counter(tuple(map(str,info(randpart(sizes)))) for _ in range(30))
    print(sizes,"dim formula 12-d-[all even] =",12-len(sizes)-(all(s%2==0 for s in sizes)),"observed (dim,weights):",dict(res))
