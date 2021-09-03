#!/usr/bin/sage

from .sboxU_cython import *

def int_to_list(i,p,n):
    l = [0 for k in range(0,n)]
    valeur = i
    for k in range(n,0,-1):
        reste = valeur % p
        l[k-1] = reste
        valeur = (valeur - reste) // p
    return l

def list_to_int(l,p,n):
    i = 0
    j = 1
    for k in range(n,0,-1):
        i += l[k-1] * j
        j *= p
    return i

def sbox_build(l,p,t,u):
    indexes = [int_to_list(i,p,t) for i in range(p**t)]
    values = [int_to_list(l[i],p,u) for i in range(p**t)]
    return PyFptFunction(p,t,u,indexes,values)

