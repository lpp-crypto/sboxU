# -*- python -*-


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

    # !TODO! implement the python class BinLinearMap
    
    def __init__(self):
        self.cpp_blm = new cpp_BinLinearMap(<std_vector[BinWord]>[ ])


    def get_input_length(self):
        return self.cpp_blm[0].get_input_length()

    
    def get_output_length(self):
        return self.cpp_blm[0].get_output_length()


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
        # !TODO!  inversion of BinLinearMap
        raise NotImplemented("still a work in progress")

    
    def rank(self):
        return self.cpp_blm[0].rank()


    def __str__(self):
        # !TODO! proper __str__ 
        return str(list(self.cpp_blm[0].get_image_vectors()))
    
    # !TODO! __rich_repr__

    # !TODO! from_blob / to_blob



def identity_BinLinearMap(int64_t n):
    return Blm([(1 << i) for i in range(0, n)])


def zero_BinLinearMap(int64_t n):
    return Blm([0 for i in range(0, n)])


def Blm(x):
    if isinstance(x, (BinLinearMap)):
        return x
    else:
        result = BinLinearMap()
        if isinstance(x, (list)):
            result.cpp_blm[0] = cpp_BinLinearMap(<std_vector[BinWord]>x)
        else:
            # !TODO! implement Blm processing of SAGE matrices 
            raise NotImplemented("Blm function cannot process input of type {}".format(type(x)))
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

            



    
