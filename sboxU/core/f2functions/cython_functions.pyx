# -*- python -*-

from sage.all import GF, Polynomial, PolynomialRing
from sage.all import Integer as SAGE_INTEGER
from sage.rings.finite_rings.finite_field_base import FiniteField

from sboxU.core.f2functions.field_arithmetic import i2f_and_f2i

from cython.operator cimport dereference



# !SECTION! Bit-fiddling

# !SUBSECTION! Wrapped C++

def oplus(BinWord x, BinWord y) -> BinWord:
    """Essentially a wrapper for the operation `^` in C++. Its purpose is to ensure that a XOR is performed regardless of the extension of the script.

    Args:
        x(BinWord): a positive integer
        y(BinWord): a positive integer

    Returns:
        A positive integer equal to the XOR of `x` and `y`.

    """
    return cpp_oplus(x, y)


def hamming_weight(BinWord x) -> int:
    """Ultimately call a C++ intrinsic to return the Hamming weight of the vector corresponding to the binary representation of `x`.
    
    Args:
        x(BinWord): a positive integer

    Returns:
        The number of bits set to 1 in the binary representation of `x`.

    """
    return cpp_hamming_weight(x)


def scal_prod(BinWord x, BinWord y) -> BinWord:
    """The canonical scalar product in F_2. Wraps a C++ function relying on specific intrinsincs.

    Args:
        x(BinWord): a positive integer
        y(BinWord): a positive integer

    Returns:
        The scalar product x⋅y, i.e. the modulo 2 sum of the products x_i y_i, where i goes from 0 to 63.
    """
    return cpp_scal_prod(x, y)


def msb(BinWord x) -> int:
    """The most significant bit.

    Args:
        x(BinWord): a positive integer

    Returns:
        The integer giving the position of the most significant bit of `x`, so that `x >> msb(x)` is always 1, unless `x` is 0. In this case, returns 0.
    """
    return cpp_msb(x)


def lsb(BinWord x) -> int:
    """The least significant bit.

    Args:
        x(BinWord): a positive integer

    Returns:
        The integer giving the position of the least significant bit set to 1 of `x`, unless `x` is 0. In this case, returns 0.
    """
    return cpp_lsb(x)

def circ_shift(BinWord x, int n, int shift) -> BinWord:
    """A circular shift is the operation of rearranging the entries in a vector, either by moving the final entry to the first position, while shifting all other entries to the next position, or by performing the inverse operation. 

    Args :
        x(BinWord) : a positive integer
        n(int) : the bit length of x 
        shift(int) : a signed integer
    Returns :
        The integer whose binary decomposition is the result of a circular shift on the binary decomposition of x by 'shift' positions. The LSB-first decomposition of x is shifted to the left if shift is positive and to the right otherwise. 
    """
    return cpp_circ_shift(x,n,shift)



# !SUBSECTION! Convenient XOR abstractions

def xor(*args) -> BinWord:
    result = 0
    for x in args:
        if isinstance(x, int):
            result = oplus(x, result)
        else:
            for y in x:
                if isinstance(y, int):
                    result = oplus(y, result)
                else:
                    raise Exception("Trying to XOR a strange type ({})".format(type(x)))
    return result



# !SECTION! Linear combinations and ranks 

 
def linear_combination(std_vector[BinWord] v, BinWord mask) -> BinWord:
    return cpp_linear_combination(v, mask)


def rank_of_vector_set(std_vector[BinWord] l) -> int:
    """Computes the rank of a set of integers interpreted as binary vectors.

    Args:
        l: a list of positive integers whose binary representation corresponds to the vector we investigate.
    Returns:
        An integer equal to the rank of the matrix obtained by concatenating these vectors. Equivalently, returns the dimension of their span.
    """
    return cpp_rank_of_vector_set(l)


# !SUBSECTION! tobin and frombin

def to_bin(BinWord x, int n) -> list:
    return cpp_to_bin(x,n)

def from_bin(std_vector[int] l) -> BinWord:
    return cpp_from_bin(l)



# !SECTION! The F2Transformation class


cdef class F2Transformation:
    """!TODO!

    Actual instanciations have to implement the following methods:
    - __getitem__(self, x: BinWord) -> BinWord that corresponds to the evaluation of the transformation on the F_2^n element corresponding to the integer x;
    - get_input_length(self) -> BinWord
    - get_output_length(self) -> BinWord
    
    """
    # !SUBSECTION! Construction and Initialization

 
    def __init__(self, name=None, input_casts : list=[], output_casts: list=[]):
        self.rename(name)
        self.input_casts = input_casts
        self.output_casts = output_casts

        
    def rename(self, name):
        if name == None:
            self.name = b"F"
        elif isinstance(name, bytes):
            self.name = name
        elif isinstance(name, str):
            self.name = name.encode("UTF-8")
        else:
            raise NotImplementedError("trying to give invalid name to S_box: {}".format(name))

        
    def attach_casts_pair(self, input_cast, output_cast) -> None:
        self.input_casts.append(input_cast)
        self.output_casts.append(output_cast)


    # !SUBSECTION! Field interaction

    def interpolate_in(self, field, var_name="X"):
        # sanity checks        
        if not isinstance(field, FiniteField):
            raise Exception("F2Transformation.interpolate_in(f) expects f to be a binary field")
        elif field.characteristic() != 2:
            raise Exception("F2Transformation.interpolate_in(f) expects f to be a binary field")
        # actual Lagrange interpolation
        i2f, f2i = i2f_and_f2i(field)
        io_pairs = [(i2f(x), i2f(self[x]))
                     for x in range(0, 2**self.get_input_length())]
        return PolynomialRing(field, var_name).lagrange_polynomial(io_pairs)
                       

    # !SUBSECTION! The __call__ method

    def __call__(self, x):
        """Querying evaluating the transformation on an input of any supported type.

        Unlike __getitem__, the input does not have to be an integer; however, it needs to be a of a type that this S_box isntance can cast to an integer. The integer obtained by querying the lookup is then cast to another type using `self.output_cast`.

        Because of the logic related to casting, it is slower than __getitem__.
        
        Args:
            x: a valid input for the cast `self.input_cast`.
        
        Returns:
            The result of calling this transformation after casting `x` to an integer, and then casting the result to the correct type.
        """
        if isinstance(x, (int, SAGE_INTEGER)):
            return self[<BinWord>x]
        else:
            for i, c in enumerate(self.input_casts):
                if c.is_valid_input(x):
                    return self.output_casts[i](self[c(x)])
            raise Exception("Could not cast input of type {} to an integer using {}".format(type(x), self.input_cast))


    
    def __getitem__(self, x):
        raise Exception("Trying to call a 'virtual' method: F2Transformation.__getitem__()")

    



