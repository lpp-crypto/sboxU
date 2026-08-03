# -*- python -*-

from sage.all import Matrix, GF, Polynomial, vector
# The following `Matrix` is not the same as above, it corresponds to an abstract type
from sage.structure.element import Matrix as SAGE_MATRIX
from sage.all import Integer as SAGE_INTEGER

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



# !SECTION! The Transformation class

cdef class F2Transformation:
    """!TODO!
    """
    # !SUBSECTION! Initialization and destruction

 
    def __init__(self, name=None, input_casts : list=[], output_casts: list=[]):
        self.rename(name)
        self.input_casts = input_casts
        self.output_casts = output_casts


        
    # !SUBSECTION! Dealing with basic attributes
    
    def rename(self, name):
        if name == None:
            self.cpp_name = new_sbox_name()
        elif isinstance(name, bytes):
            self.cpp_name = name
        elif isinstance(name, str):
            self.cpp_name = name.encode("UTF-8")
        else:
            raise NotImplementedError("trying to give invalid name to S_box: {}".format(name))

        
    def attach_casts_pair(self, input_cast, output_cast) -> None:
        self.input_casts.append(input_cast)
        self.output_casts.append(output_cast)

    

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

    def __add__(self, cst : BinWord) -> F2AffineMap:
        result = F2AffineMap()
        result.set_inner_map(dereference((<F2AffineMap>self).cpp_map) + cst)
        return result

    def __radd__(self, cst : BinWord) -> F2AffineMap:
       return self + cst

    
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
        return str(vector(GF(2),cpp_to_bin(self.get_cste(),self.get_output_length())))+ "\n+\n" + str(Matrix(
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

    def get_cste(self):
        return dereference(self.cpp_map).get_cstte()

    
    # !TODO! __rich_repr__

    # !TODO! from_blob / to_blob

    def __eq__(self,F2AffineMap L) -> bool:
        return dereference(self.cpp_map).get_image_vectors() == dereference(L.cpp_map).get_image_vectors()


# !SUBSECTION! The main factory


# !SUBSUBSECTION! Handling different subcases

def get_F2AffineMap_from_image_vectors(l, c=None, input_length=None, output_length=None) -> F2AffineMap | Exception:
    # sanity checks
    if input_length != None and len(l) != input_length:
        raise Exception("In get_F2AffineMap_from_image_vectors: mismatch between actual list length and input_length")
    if c == None:
        c = int(0)
    # ensuring compatibility
    for x in l:
        if (x >> 64) != 0:
            raise Exception("F2AffineMap image cannot be more than 64-bit long")
    # actually building the result
    result = F2AffineMap()
    result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>l, c))
    return result
    
    
def get_F2AffineMap_from_S_box(l : S_box, c=None, input_length=None, output_length=None) -> F2AffineMap | Exception:
    # sanity checks
    if input_length != None and l.get_input_length() != input_length:
        raise Exception("In get_F2AffineMap_from_S_box: mismatch between actual S_box input length and input_length")
    if output_length != None and l.get_output_length() != output_length:
        raise Exception("In get_F2AffineMap_from_S_box: mismatch between actual S_box output length and output_length")
    else:
        output_length = l.get_output_length()
    if output_length > 64:
        raise Exception("F2AffineMap image cannot be more than 64-bit long")
    if c != None and c != l[0]:
        raise Exception("In get_F2AffineMap_from_S_box: mismatch between constant and l[0]")
    # actually building the result
    result = F2AffineMap()
    result.set_inner_map(cpp_F2AffineMap(dereference((<S_box>l).cpp_sb)))
    return result

    
def get_F2AffineMap_from_univariate_Polynomial(l : Polynomial, c=None, input_length=None, output_length=None) -> F2AffineMap | Exception:
    # dealing with the underlying field 
    field = l.base_ring()
    if field.characteristic() != 2:
        raise Exception("In get_F2AffineMap_from_univariate_Polynomial: a polynomial in characteristic 2 is needed")
    n = field.degree()
    if input_length != None and input_length != n:
        raise Exception("In get_F2AffineMap_from_univariate_Polynomial: mismatch between input_size and field degree")
    if n > 64:
        raise Exception("F2AffineMap image cannot be more than 64-bit long")
    i2f, f2i = i2f_and_f2i(field)
    # dealing with the constant
    if c != None and c != f2i(l[0]):
        raise Exception("In get_F2AffineMap_from_univariate_Polynomial: mismatch between l[0] and the constant given")
    elif c == None:
        c = f2i(l[0])
    # building the mapping
    imgs = [oplus(c, f2i(l(i2f(1 << i)))) for i in range(0, n)] # we need to remove the constant part from the polynomial images
    result = F2AffineMap()
    result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>imgs, c))
    return result
    

def get_F2AffineMap_from_Matrix(l, c=None, input_length=None, output_length=None) -> F2AffineMap | Exception:
    # checking the input, and making sure that input_length, output_length and block_length are correctly set
    if isinstance(l, SAGE_MATRIX):
        # -- field coherence
        # !TODO! handle the case of a non-trivial underlying field using casts.
        # !TODO! add casts to F2AffineMap
        field = l.base_ring()
        if field.characteristic() != 2:
            raise Exception("In get_F2AffineMap_from_Matrix: the matrix must have an underlying field of characteristic 2")
        # -- lengths verifications
        block_length = field.degree()
        if input_length != None and input_length != l.ncols()*block_length:
            raise Exception("In get_F2AffineMap_from_Matrix: mismatch between l.ncols()*block_length and input_length")
        else:
            input_length = l.ncols()*block_length
        if output_length != None and output_length != l.nrows()*block_length:
            raise Exception("In get_F2AffineMap_from_Matrix: mismatch between l.nrows()*block_length and output_length")
        else:
            output_length = l.nrows()*block_length
    elif isinstance(l, (list, tuple)):
        field = GF(2)
        block_length = 1 # only GF(2) is supported in this context
        # -- checking rows
        for row in l:
            if not isinstance(row, (list, tuple)):
                raise Exception("In get_F2AffineMap_from_Matrix: the list or tuple elements must themselves be lists or tuples")
            if len(row) != len(l[0]):
                raise Exception("In get_F2AffineMap_from_Matrix: inconsistent row lengths")
        # -- checking input_length and output_length coherence
        if input_length != None and input_length != len(l[0]):
            raise Exception("In get_F2AffineMap_from_Matrix: mismatch between len(l[0]) and input_length")
        else:
            input_length = len(l[0])
        if output_length != None and output_length != len(l):
            raise Exception("In get_F2AffineMap_from_Matrix: mismatch between len(l) and output_length")
        else:
            output_length = len(l)
    # ensuring compatibility
    if output_length > 64:
        raise Exception("F2AffineMap image cannot be more than 64-bit long")
    # building image vector
    imgs = []
    if block_length == 1:
        # -- F_2 case
        for i in range(0, input_length):
            y = 0
            for j in range(0, output_length):
                if l[j][i] == 1:
                    y |= (1 << j)
            imgs.append(y)
    else:
        # -- F_{2^n} case
        i2f, f2i = i2f_and_f2i(field)
        n_blocks = input_length / block_length
        for i in range(0, input_length):
            x = i2f(1 << (i % block_length))
            y = [f2i(l[i][j] * x) for j in range(0, n_blocks)]
            y_bin = []
            for y_j in y:
                y_bin += to_bin(y_j, block_length)
            imgs.append(y_bin)
    # final steps
    result = F2AffineMap()
    result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>imgs, c))
    return result
    
    

# !SUBSUBSECTION! The main factory itself

F2AFFINEMAP_TYPE_TO_FACTORY = {
    list   : get_F2AffineMap_from_image_vectors,
    tuple  : get_F2AffineMap_from_image_vectors,
    S_box  : get_F2AffineMap_from_S_box,
    Polynomial : get_F2AffineMap_from_univariate_Polynomial,
    SAGE_MATRIX : get_F2AffineMap_from_Matrix
}

def get_F2AffineMap(l, c=0, input_length=None, output_length=None) -> F2AffineMap | Exception:

    # !TODO! rewrite this factory in the style of the get_sbox factory 
    # !TODO! add processing of matrices
    # !TODO! add processing of the offset


    if isinstance(l, (F2AffineMap)):
        return l
    
    # sanitizing because SAGE can be annoying
    if isinstance(c, SAGE_INTEGER):
        c = int(c)
    t = type(l)
    if t in F2AFFINEMAP_TYPE_TO_FACTORY.keys():
        return F2AFFINEMAP_TYPE_TO_FACTORY[t](l, c, input_length, output_length)
    elif isinstance(l, Polynomial):
        # `Polynomial` is not a real type, so indexing by it will not work
        return F2AFFINEMAP_TYPE_TO_FACTORY[Polynomial](l, c, input_length, output_length)
    elif isinstance(l, SAGE_MATRIX):
        # `SAGE_MATRIX` is not a real type, so indexing by it will not work either
        return F2AFFINEMAP_TYPE_TO_FACTORY[SAGE_MATRIX](l, c, input_length, output_length)
    else:
        raise Exception("get_F2AffineMap cannot process input of this type")




# !SUBSECTION! Common particular cases

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

def F2AffineMap_from_blocks(A, B, C, D) -> F2AffineMap:
    Ablm = get_F2AffineMap(A)
    Bblm = get_F2AffineMap(B)
    Cblm = get_F2AffineMap(C)
    Dblm = get_F2AffineMap(D)
    result = F2AffineMap()
    result.set_inner_map(cpp_F2AffineMap_from_blocks(
        dereference((<F2AffineMap>Ablm).cpp_map),
        dereference((<F2AffineMap>Bblm).cpp_map),
        dereference((<F2AffineMap>Cblm).cpp_map),
        dereference((<F2AffineMap>Dblm).cpp_map),
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

