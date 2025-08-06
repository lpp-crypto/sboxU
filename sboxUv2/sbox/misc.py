"""This module contains pure python methods to generate simple S_box
instances.

"""


from .cython_functions import S_box, Sb
from sboxUv2.f2functions import i2f_and_f2i

from random import shuffle, randint

from sage.all import GF, Integer

# !SECTION! Random SBoxes


def random_permutation_S_box(bit_length, name=None):
    """Uses the standard `shuffle` function to generate a random bijective `S_box` instance.

    Args:
        bit_length: the bit-length of the input (and output) of the function.
        name: a string intended to label the output.
    
    Returns:
        An `S_box` instance picked uniformly at random from the set of all permutations operating on the set {0, .., 2**bit_length-1}.
    
    """
    lut = list(range(0, 1 << bit_length))
    shuffle(lut)
    if name == None:
        name = "RandPerm"
    return Sb(lut, name=name)


def random_function_S_box(input_bit_length, output_bit_length, name=None):
    """Uses the standard `randint` function to generate a random `S_box` instance that is very unlikely to be bijective.

    Args:
        input_bit_length: the bit-length of the input of the function.
        output_bit_length: the bit-length of its output.
        name: a string intended to label the output.
    
    Returns:
        An `S_box` instance obtained by picking each output uniformly at random in the set {0, .., 2**output_bit_length-1}.
    
    """
    output_space_size = 1 << output_bit_length
    input_space_size  = 1 << input_bit_length
    lut = [
        randint(0, output_space_size-1)
        for x in range(0, input_space_size)
    ]
    if name == None:
        name = "RandFunc"
    return Sb(lut, name=name)



# !SECTION! Common field operations as S-boxes



def F2_mul(coeff, field=None):
    """Returns an S_box containing the lookup table of a multiplication in an extension of F_2.
    
    Args:
        coeff: the coefficient by which to multiply. Can be a field element or an integer. If an integer, then the field used must be specified.
        field: the field in which the multiplication must be made. If unspecified, the parent field of `coeff` is used.
    
    """
    if isinstance(coeff, (int, Integer)):
        if field == None:
            raise Exception("If `c` is an integer then the field must be speficied!")
        else:
            i2f, f2i = i2f_and_f2i(field)
            c = i2f(coeff)
    else:
        field = c.parent()
        i2f, f2i = i2f_and_f2i(field)
        c = coeff
    return Sb(
        [
            f2i(i2f(x) * c)
            for x in range(0, field.cardinality())
        ],
        name="Mul_{}".format(f2i(c))
    )

    

def monomial(d, field):
    """Returns an `S_box` containing the LUT of a monomial operating on the given field.

    Args:
        d: the exponent of the monomial (an integer)
        field: a finite field instance assumed to be of characteristic 2.
    """
    assert field.characteristic() == 2
    assert isinstance(d, (int, Integer))
    i2f, f2i = i2f_and_f2i(field)
    return Sb([f2i(i2f(x)**d) for x in range(0, field.cardinality())])
