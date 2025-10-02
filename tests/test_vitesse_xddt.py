from sboxUv2.statistics import xddt
from sboxUv2.core.f2functions import oplus
from time import time
from sage.crypto.sboxes import sboxes

def init_xddt(S):
    XDDT=[[[] for j in range(len(S))]for i in range(len(S))]
    for a in range(0,len(S)):
        for x in range(0,len(S)):
            XDDT[a][oplus(S[x],S[oplus(x,a)])].append(x)
    return XDDT


debut=time()
for k in sboxes.keys():
    xddt(sboxes[k])
print("Time taken with sboxUv2",time()-debut)

debut=time()
for k in sboxes.keys():
    S=list(sboxes[k])
    init_xddt(S)
print("Time taken without sboxUv2",time()-debut)
