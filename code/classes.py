# Classes of planar difference families under multipliers and translations of each block.
# usage: python classes.py P families.txt [more.txt]   (lines as printed by search: "A | B | C")
import sys
def canon(line,p):
    bl=[list(map(int,x.split())) for x in line.split('|')]
    bl.sort(key=len)
    best=None
    for u in range(1,p):
        r=[]
        for B in bl:
            m=min(tuple(sorted((u*x-u*t)%p for x in B)) for t in B); r.append(m)
        r=tuple(r)
        if best is None or r<best: best=r
    return best
def ok(line,p):
    d=sorted((x-y)%p for B in line.split('|') for x in map(int,B.split()) for y in map(int,B.split()) if x!=y); return d==list(range(1,p))
p=int(sys.argv[1]); S=[]
for f in sys.argv[2:]:
    L=[l for l in open(f) if l.strip()]; assert all(ok(l,p) for l in L); S.append({canon(l,p) for l in L}); print(f,len(L),"classes",len(S[-1]))
if len(S)>1: print("equal",S[0]==S[1])
