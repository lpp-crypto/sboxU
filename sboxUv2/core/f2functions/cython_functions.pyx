# -*- python -*-

from sage.all import Matrix, GF, Polynomial
from sboxUv2.core.f2functions.field_arithmetic import i2f_and_f2i

# !SECTION! Wrapped C++

# !SUBSECTION! Bit-fiddling

def oplus(BinWord x, BinWord y):
    """Essentially a wrapper for the operation `^` in C++. Its purpose is to ensure that a XOR is performed regardless of the extension of the script.

    Args:
        - x: a positive integer
        - y: a positive integer

    Returns:
        A positive integer equal to the XOR of `x` and `y`.

    """
    return cpp_oplus(x, y)


def hamming_weight(BinWord x):
    """Ultimately call a C++ intrinsic to return the Hamming weight of the vector corresponding to the binary representation of `x`.
    
    Args:
        - x: a positive integer

    Returns:
        The number of bits set to 1 in the binary representation of `x`.

    """
    return cpp_hamming_weight(x)


def scal_prod(BinWord x, BinWord y):
    """The canonical scalar product in F_2. Wraps a C++ function relying on specific intrinsincs.

    Args:
        - x: a positive integer
        - y: a positive integer

    Returns:
        The scalar product x⋅y, i.e. the modulo 2 sum of the products x_i y_i, where i goes from 0 to 63.
    """
    return cpp_scal_prod(x, y)


def msb(BinWord x):
    return cpp_msb(x)


def lsb(BinWord x):
    return cpp_lsb(x)



# !SUBSECTION! Linear combinations and ranks 

 
def linear_combination(std_vector[BinWord] v, BinWord mask):
    return cpp_linear_combination(v, mask)


def rank_of_vector_set(std_vector[BinWord] l):
    return cpp_rank_of_vector_set(l)


# !SUBSECTION! tobin and frombin

def to_bin(BinWord x, int n):
    return cpp_to_bin(x,n)

def from_bin(std_vector[int] l):
    return cpp_from_bin(l)

# !SUBSECTION! The BinLinearMap class

cdef class BinLinearMap:
    
    def __init__(self):
        self.cpp_blm = new cpp_BinLinearMap(<std_vector[BinWord]>[ ])


    def get_input_length(self):
        return self.cpp_blm[0].get_input_length()

    
    def get_output_length(self):
        return self.cpp_blm[0].get_output_length()

    
    def transpose(self):
        result = BinLinearMap()
        result.cpp_blm[0] = self.cpp_blm[0].transpose()
        return result

    
    def __call__(self, BinWord x):
        return self.cpp_blm[0].call(x)

    
    def __add__(self, BinLinearMap L):
        result = BinLinearMap()
        result.cpp_blm[0] = self.cpp_blm[0].add(L.cpp_blm[0])
        return result

    
    def __mul__(self, BinLinearMap L):
        result = BinLinearMap()
        result.cpp_blm[0] = self.cpp_blm[0].mul(L.cpp_blm[0])
        return result

    
    def inverse(self):
        if self.cpp_blm[0].get_input_length() > 8:
            raise NotImplemented("still a work in progress")
        else:
            result = BinLinearMap()
            result.cpp_blm[0] = self.cpp_blm[0].inverse()
            return result

    
    def rank(self):
        return self.cpp_blm[0].rank()


    def __str__(self):
        return str(Matrix(GF(2),[cpp_to_bin(x,self.get_output_length()) for x in self.cpp_blm[0].get_image_vectors()]))


    def get_S_box(self):
        result = S_box(name="L")
        (<S_box>result).cpp_sb = new cpp_S_box()
        (<S_box>result).cpp_sb[0] = self.cpp_blm[0].get_cpp_S_box()
        return result
    
    # !TODO! __rich_repr__

    # !TODO! from_blob / to_blob

    def __eq__(self,BinLinearMap L):
        return self.cpp_blm[0].get_image_vectors()==L.cpp_blm[0].get_image_vectors()


# !SUBSECTION! Factories 

def identity_BinLinearMap(int64_t n):
    return Blm([(1 << i) for i in range(0, n)])


def zero_BinLinearMap(int64_t n):
    return Blm([0 for i in range(0, n)])


def block_diagonal_BinLinearMap(A, B):
    Ablm = Blm(A)
    Bblm = Blm(B)
    result = BinLinearMap()
    result.cpp_blm[0] = cpp_block_diagonal_BinLinearMap(
        (<BinLinearMap>A).cpp_blm[0],
        (<BinLinearMap>B).cpp_blm[0],
    )
    return result
    

def Blm(l):
    if isinstance(l, (BinLinearMap)):
        return l
    else:
        result = BinLinearMap()
        if isinstance(l, (list)):
            result.cpp_blm[0] = cpp_BinLinearMap(<std_vector[BinWord]>l)
        elif isinstance(l, (S_box)):
            result.cpp_blm[0] = cpp_BinLinearMap((<S_box>l).cpp_sb[0])
        elif isinstance(l, Polynomial):
            # case of a univariate polynomial
            field = l.base_ring()
            if field.characteristic() == 2:
                n = field.degree()
                i2f, f2i = i2f_and_f2i(field)
                imgs = [f2i(l(i2f(1 << i))) for i in range(0, n)]
                result.cpp_blm[0] = cpp_BinLinearMap(<std_vector[BinWord]>imgs)
            else:
                raise Exception("A polynomial in characteristic 2 is needed for a BinLinearMap")
        else:
            # !TODO! implement Blm processing of SAGE matrices 
            raise NotImplemented("Blm function cannot process input of type {}".format(type(l)))
        return result


# !SECTION! Convenient XOR abstractions

def xor(*args):
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

            



    
