from .cython_functions import *
from ..f2functions import *

from random import shuffle, randint

from sage.all import GF

# !SECTION! Random SBoxes


def random_permutation_SBox(space_size, name=None):
    lut = list(range(0, space_size))
    shuffle(lut)
    if name == None:
        name = "RandPerm"
    return Sb(lut, name=name)


def random_function_SBox(input_space_size, output_space_size, name=None):
    lut = [
        randint(0, output_space_size-1)
        for x in range(0, input_space_size)
    ]
    if name == None:
        name = "RandFunc"
    return Sb(lut, name=name)



# !SECTION! Common field operations as S-boxes



def F2_mul(c, alphabet):
    if isinstance(alphabet, (int, Integer)):
        g = GF(2**alphabet)
    elif "degree" in dir(alphabet): # case of a field
        g = alphabet
    else:
        raise NotImplemented("'{}' is an unknown alphabet".format(type(alphabet)))
    i2f, f2i = i2f_and_f2i(alphabet)
    c_field = i2f(c)
    return Sb(
        [
            f2i(i2f(x) * c_field)
            for x in range(0, g.cardinality())
        ],
        name="Mul_{}".format(c)
    )

    

def monomial(d, field):
    i2f, f2i = i2f_and_f2i(field)
    return Sb([f2i(i2f(x)**d) for x in range(0, field.cardinality())])
