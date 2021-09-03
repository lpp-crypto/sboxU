# -*-python-*- 
# Time-stamp: <2021-08-12 16:42:46 lperrin>

from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdint cimport int64_t, uint64_t
from sboxu_cpp cimport *

BIG_SBOX_THRESHOLD = 128
DEFAULT_N_THREADS  = 2

cdef class PyFptFunction :
    cdef CppFptFunction fpt_function

    def __cinit__(self, int64_t p, int64_t t, int64_t u, vector[vector[int64_t]] indexes, vector[vector[int64_t]] values):
        self.fpt_function = CppFptFunction(p,t,u,indexes,values)

    def __call__(self, vector[int64_t] index):
        return self.fpt_function.eval(index)

    @property
    def p(self):
        return self.fpt_function.p
    @p.setter
    def p(self, p):
        self.fpt_function.p = p

    @property
    def t(self):
        return self.fpt_function.t
    @t.setter
    def t(self, t):
        self.fpt_function.t = t

    @property
    def u(self):
        return self.fpt_function.u
    @u.setter
    def u(self, u):
        self.fpt_function.u = u
    
    def fpt_ddt(self):
        return self.fpt_function.fpt_ddt()

    def fpt_differential_spectrum(self,n_threads=None):
        if n_threads == None:
            if self.p**self.t > BIG_SBOX_THRESHOLD:
                n_threads = DEFAULT_N_THREADS
            else:
                n_threads = 1
        return self.fpt_function.fpt_differential_spectrum_fast(n_threads)
