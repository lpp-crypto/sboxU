# -*- python -*-

from sage.all import Matrix, GF, Polynomial, vector
# The following `Matrix` is not the same as above, it corresponds to an abstract type
from sage.structure.element import Matrix as SAGE_MATRIX
from sage.all import Integer as SAGE_INTEGER

from sboxU.core.f2functions import i2f_and_f2i, oplus, to_bin, casts_from_field

from cython.operator cimport dereference


    

# !SECTION! The F2AffineMap class

cdef class F2AffineMap(F2Transformation):
    """This class models a linear mapping defined over F_2. It encapsulates a C++ class, `cpp_F2AffineMap`, for speed.

    While it implements methods corresponding to matrix operations, such as `transpose`, it does not rely on a matrix representation internally. Instead, it stores the vectors corresponding to the images of the canonical basis of F_2^n, and operates on these.

    Unless you are working *on* (rather than *with* `sboxU`), do not use the constructor of this class. Instead, you should rely on the `get_F2AffineMap` factory.
    
    """
    
    def __init__(self, name=None, input_casts : list=[], output_casts: list=[]):
        super().__init__(name, input_casts, output_casts)

    
    cdef set_inner_map(self, A : cpp_F2AffineMap):
        self.cpp_map = make_unique[cpp_F2AffineMap](A)

    def is_linear(self):
        return dereference(self.cpp_map).is_linear()
        
    def get_input_length(self) -> int:
        return dereference(self.cpp_map).get_input_length()

    
    def get_output_length(self) -> int:
        return dereference(self.cpp_map).get_output_length()

    
    def __getitem__(self, x : BinWord) -> BinWord:
        return dereference(self.cpp_map)(x)

    
    def __add__(self, L : F2AffineMap) -> F2AffineMap:
        result = F2AffineMap(b"(" + self.name + b"+" + L.name + b")", [], [])
        result.set_inner_map(dereference((<F2AffineMap>self).cpp_map) + dereference((<F2AffineMap>L).cpp_map))
        return result

    def __add__(self, cst : BinWord) -> F2AffineMap:
        result = F2AffineMap(b"(" + self.name + b"+" + hex(cst).encode() + b")", [], [])
        result.set_inner_map(dereference((<F2AffineMap>self).cpp_map) + cst)
        return result

    def __radd__(self, cst : BinWord) -> F2AffineMap:
       return self + cst

    
    def __hash__(self):
        # !TODO! improve the implementation of F2AffineMap.__hash__() 
        return hash(self.get_sbox())

    
    def __mul__(self, L) -> F2AffineMap | S_box:
        if isinstance(L, F2AffineMap):
            result = F2AffineMap(b"(" + self.name + b"*" + L.name + b")", [], [])
            (<F2AffineMap>result).set_inner_map(dereference((<F2AffineMap>self).cpp_map) * dereference((<F2AffineMap>L).cpp_map))
            return result
        elif isinstance(L, S_box):
            result = S_box(b"(" + self.name + b"*" + L.name + b")", [], [])
            (<S_box>result).set_inner_sbox(
                (<cpp_S_box>dereference((<F2AffineMap>self).cpp_map).get_cpp_S_box()).mul(
                    <cpp_S_box>dereference((<S_box>L).cpp_sb)
                )
            )
            return result
        else:
            raise NotImplemented()

    
    def inverse(self) -> F2AffineMap | Exception:
        if (dereference(self.cpp_map).is_invertible()):
            result = F2AffineMap(self.name + b"^-1", self.input_casts, self.output_casts)
            result.set_inner_map(dereference(self.cpp_map).inverse())
            return result
        else:
            print("Trying to invert a non-invertible F2AffineMap")
            raise Exception("Trying to invert a non-invertible F2AffineMap")

    
    def transpose(self) -> F2AffineMap:
        result = F2AffineMap(self.name + b"^T", self.input_casts, self.output_casts)
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


    def get_sbox(self) -> S_box:
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


# !SECTION! Building F2AffineMap:s


# !SUBSECTION! Handling different subcases

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
    result = F2AffineMap(None, [], [])
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
    result = F2AffineMap(None, [], []) # !TODO! proper handling of S-box name when turning it into an F2AffineMap
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
    result = F2AffineMap(None, [], [])
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
    result = F2AffineMap(None, [], [])
    result.set_inner_map(cpp_F2AffineMap(<std_vector[BinWord]>imgs, c))
    return result
    
    

# !SUBSECTION! The main factory itself

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


