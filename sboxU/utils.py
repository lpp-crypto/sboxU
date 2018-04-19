#!/usr/bin/sage

from sage.all import *
from sboxu_cpp import oplus_cpp

def oplus(x, y):
    return oplus_cpp(long(x), long(y))
