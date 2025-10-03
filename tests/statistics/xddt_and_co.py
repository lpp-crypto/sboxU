from sboxUv2.statistics import xddt,yddt,zddt
from sage.crypto.sboxes import sboxes
Sb0=sboxes["Midori_Sb0"]

XDDT=xddt(Sb0)
YDDT=yddt(Sb0)
ZDDT=zddt(Sb0)

print("Testing the XDDT of Midori_Sb0")

if set(XDDT[0xa][0xa])=={3,4,9,14} and set(XDDT[0xa][0xf])=={1,6,11,12} and set(XDDT[0xf][0xa])=={0,5,10,15} and set(XDDT[0xf][0xf])=={2,7,8,13}:
    print("Success")
else :
    print("Failure")

print("Testing the YDDT of Midori_Sb0")

if set(YDDT[0xa][0xa])=={3,4,9,14} and set(YDDT[0xa][0xf])=={0,5,10,15} and set(YDDT[0xf][0xa])=={1,6,11,12} and set(YDDT[0xf][0xf])=={2,7,8,13}:
    print("Success")
else :
    print("Failure")

print("Testing the ZDDT of Midori_Sb0")

if set(ZDDT[0xa][0xa])=={51,78,153,228} and set(ZDDT[0xa][0xf])=={12,91,161,246} and set(ZDDT[0xf][0xa])=={26,111,181,192} and set(ZDDT[0xf][0xf])=={45,119,136,210}:
    print("Success")
else :
    print("Failure")