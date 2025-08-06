# -*- python -*-


from ..config import N_THREADS
from ..display import stylize

from sboxUv2.sbox import Sb # don't replace this with *: it is crucial to S_box not be imported in this import call
from sboxUv2.sbox.cython_functions cimport *


# !SECTION! Dealing with Spectra

# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    def __init__(self):
        self.cpp_sp = new cpp_Spectrum()


    cdef set_inner_sp(Spectrum self, cpp_Spectrum sp):
        self.cpp_sp[0] = sp
        

    def __getitem__(self, k):
        return self.cpp_sp.brackets(k)

    
    def __len__(self):
        return self.cpp_sp.size()


    def __str__(self):
        result = "{"
        ks = self.keys()
        for k in ks[:-1]:
            result += "{:d}:{:d}, ".format(k, self.cpp_sp.brackets(k))
        result += stylize(
            "{:d}:{:d}".format(
                ks[-1],
                self.cpp_sp.brackets(ks[-1])
                ),
            "b"
        )
        return result + "}"
    
    
    def keys(self):
        return self.cpp_sp.keys()

    
    def incr(self, x):
        self.cpp_sp.incr(x)


    def incr_by_counting(self, v):
        self.cpp_sp.incr_by_counting(v)

        
    def maximum(self):
        return self.cpp_sp.maximum()
    

    def occurrence_of_maximum(self):
        return self.cpp_sp.brackets(self.cpp_sp_maximum())


# !SUBSECTION! A Spectrum factory

def Spctr(x):
    if isinstance(x, list):
        result = Spectrum()
        result.incr_by_counting(x)
        return result
    else:
        return NotImplemented("{} is not a valid input type for Spctr".format(type(x)))



# !SECTION! Differential Properties

def differential_spectrum(_s):
    s = Sb(_s)
    result = Spectrum()
    result.set_inner_sp(
        cpp_differential_spectrum((<S_box>s).cpp_sb[0],
                                  N_THREADS)
    )
    return result
    
def ddt(_s):
    s = Sb(_s)
    result = cpp_ddt((<S_box>s).cpp_sb[0])
    return result
