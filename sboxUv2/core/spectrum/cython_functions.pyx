# -*- python -*-

from sboxUv2.config import n_threads_from_sbox_size


# !SECTION! Dealing with Spectra

# !SUBSECTION! The Spectrum class

cdef class Spectrum:
    def __init__(self, name=b"Spc"):
        self.name = name
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
        for k in ks:
            result += "{:d}:{:d}, ".format(k, self.cpp_sp.brackets(k))
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
                    self.cpp_sp.brackets(k)
                )
            result += "[bold blue]{:d}[/][black]:{:d}[/], ".format(
                    ks[-1],
                    self.cpp_sp.brackets(ks[-1])
                )
            return result[:-2] + "}"

    
    def __iter__(self):
        for x in self.cpp_sp.keys():
            yield x


    def absolute(self):
        result = Spectrum(name=b"|" + self.name + b"|")
        for k in self.keys():
            result.incr_by_amount(abs(k), self.cpp_sp.brackets(k))
        return result
        
    
    def keys(self):
        return self.cpp_sp.keys()

    
    def incr(self, x):
        self.cpp_sp.incr(x)

    
    def incr_by_amount(self, x, c):
        self.cpp_sp.incr_by_amount(x, c)


    def incr_by_counting(self, v):
        self.cpp_sp.incr_by_counting(v)

        
    def maximum(self):
        return self.cpp_sp.maximum()
    

    def occurrence_of_maximum(self):
        return self.cpp_sp.brackets(self.cpp_sp_maximum())


# !SUBSECTION! A Spectrum factory


def Spctr(x):
    if isinstance(x, (list)):
        result = Spectrum()
        result.incr_by_counting(x)
        return result
    else:
        return NotImplemented("{} is not a valid input type for Spctr".format(type(x)))

