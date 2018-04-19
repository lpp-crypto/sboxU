#!/usr/bin/sage

from sboxu_cpp import oplus_cpp

def oplus(x, y):
    return oplus_cpp(long(x), long(y))
