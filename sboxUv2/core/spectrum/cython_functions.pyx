# -*- python -*-

from sboxUv2.config import n_threads_from_sbox_size

from libcpp.memory cimport unique_ptr, make_unique
from cython.operator cimport dereference


# !SECTION! Dealing with Spectra

# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    def __init__(self, name=b"Spc"):
        self.name = name
        self.cpp_sp = make_unique[cpp_Spectrum]()
       
        

    cdef set_inner_sp(Spectrum self, cpp_Spectrum sp):
        self.cpp_sp = make_unique[cpp_Spectrum](sp)
        

    def __getitem__(self, int64_t k):
        return dereference(self.cpp_sp)[k]

    
    def __len__(self):
        return dereference(self.cpp_sp).size()


    def __str__(self):
        result = "{"
        ks = self.keys()
        for k in ks:
            result += "{:d}:{:d}, ".format(k, dereference(self.cpp_sp)[k])
        return result[:-2] + "}"


    def __rich_str__(self):
        result = "[bright_black]{:12s}[/] ".format(self.name.decode("UTF-8"))
        result += "{"
        ks = self.keys()
        if len(ks) == 0:
            result += "}"
            return result
        else:
            for k in ks[:-1]:
                result += "[bold bright_black]{:d}[/][black]:[/][white]{:d}[/], ".format(
                    k,
                    dereference(self.cpp_sp)[k]
                )
            result += "[bold blue]{:d}[/][black]:{:d}[/], ".format(
                    ks[-1],
                    dereference(self.cpp_sp)[ks[-1]]
                )
            return result[:-2] + "}"

    
    def __iter__(self):
        for x in dereference(self.cpp_sp).cpp_keys():
            yield x


    def absolute(self):
        result = Spectrum(name=b"|" + self.name + b"|")
        result.set_inner_sp(dereference(self.cpp_sp).absolute())
        return result
        
    
    def keys(self):
        return dereference(self.cpp_sp).cpp_keys()

    
    def incr(self, x):
        dereference(self.cpp_sp).incr(x)

    
    def incr_by_amount(self, x, c):
        dereference(self.cpp_sp).incr_by_amount(x, c)


    def incr_by_counting(self, v):
        dereference(self.cpp_sp).incr_by_counting(v)

        
    def maximum(self):
        return dereference(self.cpp_sp).maximum()
    

    def occurrence_of_maximum(self):
        return dereference(self.cpp_sp)[self.cpp_sp_maximum()]


# !SUBSECTION! A Spectrum factory


def get_Spectrum(x):
    if isinstance(x, (list)):
        result = Spectrum()
        result.incr_by_counting(x)
        return result
    else:
        return NotImplemented("{} is not a valid input type for get_Spectrum".format(type(x)))

