# -*- python -*-



from sboxUv2.f2functions import *
from sboxUv2.f2functions.cython_functions cimport *

from sage.crypto.sboxes import SBox as sage_SBox
from sage.all import Integer


# !SECTION! Helpers


# !SUBSECTION! Functions

cdef cpp_S_box pyx_add_sboxes(cpp_S_box s, cpp_S_box t):
    """A wrapper for the cpp_S_box overloaded operator +."""
    return s.add(t)

cdef cpp_S_box pyx_mul_sboxes(cpp_S_box s, cpp_S_box t):
    """A wrapper for the cpp_S_box overloaded operator *."""
    return s.mul(t)


def new_sbox_name():
    """Returns a unique name that can be given to an S-box.

    It uses the module variable `sboxU_SBOXES_COUNTER` to achieve this goal by incrementing it each time it is used.

    Returns:
        A bytearray corresponding to the next unique S_box name.
    """
    global sboxU_SBOXES_COUNTER
    result = bytes("S_{}".format(sboxU_SBOXES_COUNTER).encode("UTF-8"))
    sboxU_SBOXES_COUNTER += 1
    return result


    
# !SUBSECTION! Global variables 

# global variable used to give unique names to S-boxes
cdef uint64_t sboxU_SBOXES_COUNTER = 0


# !SECTION! The S_box class



cdef class S_box:
    # "cdef" attributes are declared in the .pxd file
    """The S_box class stores the lookup table of an vectorial boolean function, and provides useful methods to interact with it.

    Objects of this class should be initialized using the :py:func:Sb function.
    """
                                 
    
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
            raise NotImplemented("trying to give invalid name to S_box: {}".format(name))

        
    # !SUBSECTION! Python built-in methods
    
    def __add__(S_box self, _s):
        """Addition in F_2 (i.e., XOR).

        Args:
            _s: the S_box to add to the current one. Must be an S_boxable type.
        
        Returns:
            An `S_box` instance whose output is the XOR of `self` and `_s`. 
        """
        s = Sb(_s)
        if len(s) != len(self):
            raise Exception("Trying to add S_boxes of different lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the S_box.cpp_sb pointer
        name = self.cpp_name + b"+" + s.name()
        result = S_box(name)
        (<S_box>result).set_inner_sbox(pyx_add_sboxes(self.cpp_sb[0], (<S_box>s).cpp_sb[0]))
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
        """Querying the S-box on a specific input.
        
        Args:
            x: an integer whose binary representation corresponds to the input on which to query the S-box.
        
        Returns:
            The result of calling this S-box on the input of `x`.
        """
        return self.cpp_sb.brackets(x)


    def __len__(self):
        """Returns:
            The number of entries in the lookup table of this S_box.
        """
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
        """Composing the S_box with itself (or with its inverse).

        Args:
            d: the number of times that the function should be iterated.

        Returns:
            If `d` is equal to 0, returns the identity S_Box. If it is a positive integer, returns the S_box corresponding to d iterations of the current S_box. If it is negative, does the same but for its inverse (throws an exception if the S_box is not invertible).
        
        The `modulo` argument is needed by the python syntax for the __pow__ function. An error will be thrown if it is set.
        """
        if modulo != None:
            raise NotImplemented("why are you using a modulo (second pow argument) here?")
        elif d == 0:
            return identity_S_box(len(self))
        elif d == 1:
            return self
        elif d == -1:
            return self.inverse()
        else:
            result = S_box(name=self.name() + b"**" + str(d).encode("UTF-8"))
            (<S_box>result).set_inner_sbox(cpp_S_box(<cpp_vector[uint64_t]> list(self.input_space())))
            if d >= 0:
                for i in range(0, d):
                    (<S_box>result).cpp_sb[0] = pyx_mul_sboxes((<S_box>self).cpp_sb[0],
                                                              (<S_box>result).cpp_sb[0])
                return result
            elif self.is_invertible():
                inv = self.inverse()
                for i in range(0, -d):
                    (<S_box>result).cpp_sb[0] = pyx_mul_sboxes((<S_box>inv).cpp_sb[0],
                                                              (<S_box>result).cpp_sb[0])
                return result
            else:
                raise Exception("Trying to compute the negative power of a non-bijective function")


    def __lmul__(self):
        return NotImplemented("Only right composition is implemented for objects of class S_box")

    
    def __mul__(S_box self, _s):
        """The composition operator.

        Args:
            _s: an S_boxable object with which the current S_box must be right-composed.

        Returns:
            An S_box corresponding to the function "F o _s", where F is the current S-box, and _s is the input to the function.
        """
        s = Sb(_s)
        if self.get_input_length() != s.get_output_length():
            raise Exception("Trying to compose S_boxes of incompatible lengths:\n{}\n{}".format(self, s))
        # below, the [0] is used to follow the S_box.cpp_sb pointer
        name = self.cpp_name + "◦".encode("UTF-8") + s.name()
        result = S_box(name)
        (<S_box>result).set_inner_sbox(pyx_mul_sboxes(self.cpp_sb[0], (<S_box>s).cpp_sb[0]))
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


    cdef set_inner_sbox(S_box self, cpp_S_box s):
        if self.cpp_sb:
            del self.cpp_sb
        self.cpp_sb = new cpp_S_box()
        self.cpp_sb[0] = s


    def name(self):
        return self.cpp_name

    # !SUBSUBSECTION! Components and coordinates
    
    def coordinate(S_box self, uint64_t i):
        """Args:
            i: the index of the coordinate, where 0 is the bit of lowest weight.
        
        Returns:
            An S_box instance mapping n bits to 1 corresponding to the i-th coordinate of S.
        
        """
        assert i < self.cpp_sb.get_output_length()
        result = S_box(name=self.cpp_name + ("_{:x}".format(i)).encode("UTF-8"))
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.coordinate(<uint64_t>i))
        return result
        
    
    def component(S_box self, uint64_t a):
        """Returns:
            An S_box instance mapping n bits to 1 corresponding to the component x \mapsto S(x) \cdot a, where \cdot is the standard scalar product.
        
        """
        result = S_box(name=self.cpp_name + ("•{:x}".format(a)).encode("UTF-8"))
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.component(<uint64_t>a))
        return result
        

    # !SUBSUBSECTION! Derivatives

    def derivative(S_box self, uint64_t delta):
        """Returns:
            An S_box of the same dimension as S corresponding to its derivative in the direction delta, i.e. x \mapsto S(x+delta)+S(x).
        
        """
        result = S_box(name=("Δ_{:x} ".format(delta)).encode("UTF-8") + self.cpp_name)
        result.set_inner_sbox(<cpp_S_box>self.cpp_sb.derivative(<uint64_t>delta))
        return result
        
    
    
    # !SUBSECTION! Function composition

    def is_invertible(self):
        """Returns:
            True if the current S_box is a bijection, False otherwise.
        """
        return self.cpp_sb.is_invertible()

    
    def inverse(self):
        """Returns:
            An S_box instance corresponding to the compositional inverse of the current S_box.

        If the current S_box is not invertible, will probably crash.
        """
        name = self.cpp_name + b"^-1"
        result = S_box(name=name)
        (<S_box>result).set_inner_sbox(<cpp_S_box>(self.cpp_sb.inverse()))
        return result



# !SECTION! Generating S_boxes


# !SUBSECTION! Main factory 

def Sb(s, name=None):
    """Turns its input into an object of the S_box class.

    If it is already an S_box instance, simply returns its
    input. Otherwise, builds the lookup table, and then create the
    corresponding S_box instance.

    Args:
        - s: an object of a class that can be turned into an S_box.
        - name: the name to give the object. If none is provided, one
          will be picked using `sboxU_SBOXES_COUNTER`.

    """
    if isinstance(s, S_box):
        return s
    else:
        result = S_box(name=name)
        if isinstance(s, list):
            (<S_box>result).cpp_sb = new cpp_S_box(<cpp_vector[uint64_t]>s)
        elif isinstance(s, sage_SBox):
            (<S_box>result).cpp_sb = new cpp_S_box(<cpp_vector[uint64_t]>list(s))
            # !TODO! add other possible initializations (from polynomials for ex) 
        else:
            try:
                result = s.get_S_box()
            except:
                raise NotImplemented("can't turn object of type '{}' into an S_box".format(type(s)))
    return result


# !SUBSECTION! Simple structures

cdef S_box pyx_F2_trans(uint64_t k, n):
    """Wrapper for the `cpp_translation` function. """
    result = S_box(name="Add_{}".format(k))
    result.cpp_sb = new cpp_S_box()
    result.cpp_sb[0] = cpp_translation(k, n)
    return result


def F2_trans(additive_cstte, field=None, bit_length=None):
    """Returns an S_box containing the lookup table of a simple XOR over a given field extension of F_2.

    If additive_cstte is an integer, then either `field` or `bit_length` must be set. If it is a field element, both `field` and `bit_length` will be ignored.
    
    Args:
        additive_cstte: the constant to add. Can be a field element or an integer. If an integer, then the field used must be specified.
        field: the field in which the multiplication must be made if `additive_cstte` is an integer.
        bit_length: the bit-length to use for both the input and output if `additive_cstte` is an integer.

    Returns:
        An S_box instance
    """
    if isinstance(additive_cstte, (int, Integer)):
        k = additive_cstte
        if isinstance(bit_length, (int, Integer)):
            n = bit_length
        elif "degree" in dir(field): # case of a field
            n = field.degree()
        else:
            raise Exception("If `additive_cstte` is an integer then either `field` or `bit_length` must be speficied!")
    else: # case where the additive constant is a finite field element
        k = ffe_to_int(additive_cstte)
        n = additive_cstte.parent().degree()
    return pyx_F2_trans(k, n)


def identity_S_box(length):
    """Returns an S_box instance corresponding to the identity
    function, i.e. the one mapping x to itself.

    """
    return Sb(list(range(0, length)))
