# -*- python -*-

from ..f2functions import *
from ..f2functions.cython_functions cimport *

from sage.crypto.sboxes import SBox as sage_SBox


# !SECTION! Helpers


# !SUBSECTION! Functions

cdef cpp_SBox pyx_add_sboxes(cpp_SBox s, cpp_SBox t):
    return s.add(t)

cdef cpp_SBox pyx_mul_sboxes(cpp_SBox s, cpp_SBox t):
    return s.mul(t)


def new_sbox_name():
    global sboxU_SBOXES_COUNTER
    result = bytes("S_{}".format(sboxU_SBOXES_COUNTER).encode("UTF-8"))
    sboxU_SBOXES_COUNTER += 1
    return result


    
# !SUBSECTION! Global variables 

# global variable used to give unique names to S-boxes
cdef uint64_t sboxU_SBOXES_COUNTER = 0


# !SECTION! The SBox class



cdef class SBox:
    # "cdef" attributes are declared in the .pxd file
                                 
    
    # !SUBSECTION! Initialization

 
    def __init__(self, name=None):
        self.rename(name)


    # !SUBSECTION! Dealing with the name
    
    def rename(self, name):
        if name == None:
            self.cpp_name = new_sbox_name()
        elif isinstance(name, bytes):
            self.cpp_name = name
        elif isinstance(name, str):
            self.cpp_name = name.encode("UTF-8")
        else:
            raise NotImplemented("trying to give invalid name to SBox: {}".format(name))

        
    # !SUBSECTION! Python built-in methods
    
    def __add__(SBox self, _s):
        s = Sb(_s)
        if len(s) != len(self):
            raise Exception("Trying to add SBoxes of different lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the SBox.cpp_sb pointer
        name = self.cpp_name + b"+" + s.name()
        result = SBox(name)
        (<SBox>result).set_inner_sbox(pyx_add_sboxes(self.cpp_sb[0], (<SBox>s).cpp_sb[0]))
        return result
        

    def __eq__(self, s):
        if len(s) != self.cpp_sb.size():
            return False
        else:
            for x in self.input_space():
                if s[x] != self[x]:
                    return False
        return True


    def __ne__(self, s):
        return not self.__equal__(s)

        
    def __getitem__(self, uint64_t x):
        return self.cpp_sb.brackets(x)


    def __len__(self):
        return self.cpp_sb.size()

    
    def __str__(self):
        return "({:2d},{:2d}) {} = {}".format(
            self.get_input_length(),
            self.get_output_length(),
            self.cpp_name.decode("UTF-8"),
            self.cpp_sb.content_string_repr().decode("UTF-8")
        )

    
    def __iter__(self):
        for x in self.input_space():
            yield self[x]

    
    def __hash__(self):
        return hash(self.cpp_sb.content_string_repr())



    def __pow__(self, d, modulo):
        if modulo != None:
            raise NotImplemented("why are you using a modulo (second pow argument) here?")
        elif d == 0:
            return identity_SBox(len(self))
        elif d == 1:
            return self
        elif d == -1:
            return self.inverse()
        else:
            result = SBox(name=self.name() + b"**" + str(d).encode("UTF-8"))
            (<SBox>result).set_inner_sbox(cpp_SBox(<cpp_vector[uint64_t]> list(self.input_space())))
            if d >= 0:
                for i in range(0, d):
                    (<SBox>result).cpp_sb[0] = pyx_mul_sboxes((<SBox>self).cpp_sb[0],
                                                              (<SBox>result).cpp_sb[0])
                return result
            elif self.is_invertible():
                inv = self.inverse()
                for i in range(0, -d):
                    (<SBox>result).cpp_sb[0] = pyx_mul_sboxes((<SBox>inv).cpp_sb[0],
                                                              (<SBox>result).cpp_sb[0])
                return result
            else:
                raise Exception("Trying to compute the negative power of a non-bijective function")


    def __lmul__(self):
        return NotImplemented("Only right composition is implemented for objects of class SBox")

    
    def __mul__(SBox self, _s):
        s = Sb(_s)
        if self.get_input_length() != s.get_output_length():
            raise Exception("Trying to compose SBoxes of incompatible lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the SBox.cpp_sb pointer
        name = self.cpp_name + b"*" + s.name()
        result = SBox(name)
        (<SBox>result).set_inner_sbox(pyx_mul_sboxes(self.cpp_sb[0], (<SBox>s).cpp_sb[0]))
        return result
            
        
    # !SUBSECTION! Functions dealing with input/output sizes

    
    def get_input_length(self):
        return self.cpp_sb.get_input_length()

    
    def input_space_size(self):
        return self.cpp_sb.input_space_size()


    def input_space(self):
        return range(0, self.cpp_sb.input_space_size())

    
    def get_output_length(self):
        return self.cpp_sb.get_output_length()


    def output_space_size(self):
        return self.cpp_sb.output_space_size()

    
    def output_space(self):
        return range(0, self.cpp_sb.output_space_size())

    
    # !SUBSECTION! Basic accessors
    
    def lut(self):
        return self.cpp_sb.get_lut()


    cdef set_inner_sbox(SBox self, cpp_SBox s):
        if self.cpp_sb:
            del self.cpp_sb
        self.cpp_sb = new cpp_SBox()
        self.cpp_sb[0] = s


    def name(self):
        return self.cpp_name
    
    
    # !SUBSECTION! Function composition

    def is_invertible(self):
        return self.cpp_sb.is_invertible()

    
    def inverse(self):
        name = self.cpp_name + b"^-1"
        result = SBox(name=name)
        (<SBox>result).set_inner_sbox(<cpp_SBox>(self.cpp_sb.inverse()))
        return result



# !SECTION! Generating SBoxes


# !SUBSECTION! Main factory 

def Sb(s, name=None):
    """Turns its input into an object of the SBox class.

    If it is already an SBox instance, simply returns its
    input. Otherwise, builds the lookup table, and then create the
    corresponding SBox instance.

    Args:
        - s: an object of a class that can be turned into an SBox.
        - name: the name to give the object. If none is provided, one
          will be picked using `sboxU_SBOXES_COUNTER`.

    """
    if isinstance(s, SBox):
        return s
    else:
        result = SBox(name=name)
        if isinstance(s, list):
            (<SBox>result).cpp_sb = new cpp_SBox(<cpp_vector[uint64_t]>s)
        elif isinstance(s, sage_SBox):
            (<SBox>result).cpp_sb = new cpp_SBox(<cpp_vector[uint64_t]>list(s))
            # !TODO! add other possible initializations (from polynomials for ex) 
        else:
            try:
                result = s.get_SBox()
            except:
                raise NotImplemented("can't turn object of type '{}' into an SBox".format(type(s)))
    return result


# !SUBSECTION! Simple structures

cdef SBox pyx_F2_trans(uint64_t k, alphabet):
    """Returns a translation by k over a finite field of
    characteristic, i.e. a function that XORs `k` into its input.

    """
    if isinstance(alphabet, (int, Integer)):
        n = alphabet
    elif "degree" in dir(alphabet): # case of a field
        n = alphabet.degree()
    else:
        raise NotImplemented("'{}' is an unknown alphabet".format(type(alphabet)))
    result = SBox(name="Add_{}".format(k))
    result.cpp_sb = new cpp_SBox()
    result.cpp_sb[0] = cpp_translation(k, n)
    return result


def F2_trans(k, alphabet):
    return pyx_F2_trans(k, alphabet)


def identity_SBox(length):
    return Sb(list(range(0, length)))
