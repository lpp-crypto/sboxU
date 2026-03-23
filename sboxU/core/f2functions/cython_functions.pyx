# -*- python -*-

from sage.all import Matrix, GF, Polynomial
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

# !SECTION! The F2AffineMap class

cdef class F2AffineMap:
    """This class models a linear mapping defined over F_2. It encapsulates a C++ class, `cpp_F2AffineMap`, for speed.

    While it implements methods corresponding to matrix operations, such as `transpose`, it does not rely on a matrix representation internally. Instead, it stores the vectors corresponding to the images of the canonical basis of F_2^n, and operates on these.

    Unless you are working *on* (rather than *with* `sboxU`), do not use the constructor of this class. Instead, you should rely on the `get_F2AffineMap` factory.
    
    """
    
    def __init__(self):
        pass

    
    cdef set_inner_map(self, A : cpp_F2AffineMap):
        self.cpp_map = make_unique[cpp_F2AffineMap](A)

    def is_linear(self):
        return dereference(self.cpp_map).is_linear()
        
    def get_input_length(self) -> int:
        return dereference(self.cpp_map).get_input_length()

    
    def get_output_length(self) -> int:
        return dereference(self.cpp_map).get_output_length()

    
    def __call__(self, x : BinWord) -> BinWord:
        return dereference(self.cpp_map)(x)

    
    def __add__(self, L : F2AffineMap) -> F2AffineMap:
        result = F2AffineMap()
        result.set_inner_map(dereference((<F2AffineMap>self).cpp_map) + dereference((<F2AffineMap>L).cpp_map))
        return result

    
    def __hash__(self):
        # !TODO! improve the implementation of F2AffineMap.__hash__() 
        return hash(self.get_S_box())

    
    def __mul__(self, F2AffineMap L) -> F2AffineMap:
        result = F2AffineMap()
        result.set_inner_map(dereference((<F2AffineMap>self).cpp_map) * dereference((<F2AffineMap>L).cpp_map))
        return result

    
    def inverse(self) -> F2AffineMap | Exception:
        if (dereference(self.cpp_map).is_invertible()):
            result = F2AffineMap()
            result.set_inner_map(dereference(self.cpp_map).inverse())
            return result
        else:
            print("Trying to invert a non-invertible F2AffineMap")
            raise Exception("Trying to invert a non-invertible F2AffineMap")

    
    def transpose(self) -> F2AffineMap:
        result = F2AffineMap()
        result.set_inner_map(dereference(self.cpp_map).transpose())
        return result

    
    def rank(self) -> int:
        return dereference(self.cpp_map).rank()


    def __str__(self) -> str:
        return str(Matrix(
            GF(2),
            [cpp_to_bin(x, self.get_output_length())
             for x in dereference(self.cpp_map).get_image_vectors()]
        ).transpose())


    def get_S_box(self) -> S_box:
        result = S_box(name="L")
        result.set_inner_sbox(dereference(self.cpp_map).get_cpp_S_box())
        return result
    
    
    def get_image_vectors(self):
        return dereference(self.cpp_map).get_image_vectors()

    
    # !TODO! __rich_repr__

    # !TODO! from_blob / to_blob

    def __eq__(self,F2AffineMap L) -> bool:
        return dereference(self.cpp_map).get_image_vectors() == dereference(L.cpp_map).get_image_vectors()


# !SUBSECTION! Factories 

def identity_F2AffineMap(int64_t n) -> F2AffineMap:
    return get_F2AffineMap([(1 << i) for i in range(0, n)],n,n)


def zero_F2AffineMap(n : BinWord, m : BinWord) -> F2AffineMap:
    return get_F2AffineMap([0 for i in range(0, n)], n, m)


def block_diagonal_F2AffineMap(A, B) -> F2AffineMap:
    Ablm = get_F2AffineMap(A)
    Bblm = get_F2AffineMap(B)
    result = F2AffineMap()
    result.set_inner_map(cpp_block_diagonal_F2AffineMap(
        dereference((<F2AffineMap>Ablm).cpp_map),
        dereference((<F2AffineMap>Bblm).cpp_map),
    ))
    return result


def circ_shift_F2AffineMap(int n, int shift) -> F2AffineMap:
    """A circular shift is the operation of rearranging the entries in a vector, either by moving the final entry to the first position, while shifting all other entries to the next position, or by performing the inverse operation. 

    Args : 
        - n : a positive integer
        - shift : a signed integer
    Returns :
        A F2AffineMap object which encodes the circular shift by 'shift' positions. This linear map is an automorphism of (F_2)^n. As for circ_shift, the LSB-first decomposition of a vector x is shifted to the left if shift is positive and to the right otherwise. 
    """
    return get_F2AffineMap([circ_shift(1 <<i,n,shift) for i in range(0, n)], n, n)

def bit_permutation_F2AffineMap(p) -> F2AffineMap:
    """
    A bit permutation is the operation of rearranging the entries in a F2 vector, according to a given permutation. 
    
    Args :
        p : The lut of the bit permutation.
    
    Returns : 
        A BinLinearMap corresponding to bit permutation associated to p.
    """
    return get_F2AffineMap([1 << p[i] for i in range(len(p))])
    
# !SUBSECTION! The main factory


def get_F2AffineMap(l, input_length=None, output_length=None) -> F2AffineMap:
    if isinstance(l, (F2AffineMap)):
        return l
    elif input_length is None and output_length is None:
        result = F2AffineMap()
        if isinstance(l, (list)):
            result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>l))
        elif isinstance(l, (S_box)):
            result.set_inner_map(cpp_F2AffineMap(dereference((<S_box>l).cpp_sb)))
        elif isinstance(l, Polynomial):
            # case of a univariate polynomial
            field = l.base_ring()
            if field.characteristic() == 2:
                n = field.degree()
                i2f, f2i = i2f_and_f2i(field)
                imgs = [f2i(l(i2f(1 << i))) for i in range(0, n)]
                result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>imgs))
            else:
                raise Exception("A polynomial in characteristic 2 is needed for a F2AffineMap")
        else:
            # !TODO! implement Blm processing of SAGE matrices 
            raise NotImplemented("Blm function cannot process input of type {}".format(type(l)))
        return result
    elif (input_length is None) or (output_length is None):
        raise NotImplemented("You must specify both input_length and output_length or none of them")
    else :
        result = F2AffineMap()
        if isinstance(l, (list)):
            result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>l,
                                                 input_length,
                                                 output_length))
        elif isinstance(l, (S_box)):
            if l.get_input_length()!=input_length or l.get_output_length()!= output_length:
                raise Exception("You specified uncompatible input or output length")
            result.set_inner_map(cpp_F2AffineMap(dereference((<S_box>l).cpp_sb)))
        elif isinstance(l, Polynomial):
            # case of a univariate polynomial
            field = l.base_ring()
            if field.characteristic() == 2:
                n = field.degree()
                i2f, f2i = i2f_and_f2i(field)
                imgs = [f2i(l(i2f(1 << i))) for i in range(0, n)]
                result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>imgs,
                                                     input_length,
                                                     output_length))
            else:
                raise Exception("A polynomial in characteristic 2 is needed for a F2AffineMap")
        else:
            # !TODO! implement Blm processing of SAGE matrices 
            raise NotImplemented("Blm function cannot process input of type {}".format(type(l)))
        return result

