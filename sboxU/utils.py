#!/usr/bin/sage

from sage.all import *
from sboxu_cpp import oplus_cpp

def oplus(x, y):
    return oplus_cpp(long(x), long(y))

def inverse(s):
    result = [0 for i in xrange(0, len(s))]
    for x in xrange(0, len(s)):
        result[s[x]] = x
    return result
    
