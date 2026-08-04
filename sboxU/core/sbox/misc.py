"""This module contains pure python methods to generate simple S_box
instances.

"""


from sboxU.core.sbox.cython_functions import S_box, get_sbox
from sboxU.core.f2functions import i2f_and_f2i

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


    # !TODO! this function shouldn't be here 
    
    """
    lut = list(range(0, 1 << bit_length))
    shuffle(lut)
    if name == None:
        name = "RandPerm"
    return get_sbox(lut, name=name)


def random_function_S_box(input_bit_length, output_bit_length, name=None):
    """Uses the standard `randint` function to generate a random `S_box` instance that is very unlikely to be bijective.

    Args:
        input_bit_length: the bit-length of the input of the function.
        output_bit_length: the bit-length of its output.
        name: a string intended to label the output.
    
    Returns:
        An `S_box` instance obtained by picking each output uniformly at random in the set {0, .., 2**output_bit_length-1}.

    # !TODO! this function shouldn't be here
    
    """
    output_space_size = 1 << output_bit_length
    input_space_size  = 1 << input_bit_length
    lut = [
        randint(0, output_space_size-1)
        for x in range(0, input_space_size)
    ]
    if name == None:
        name = "RandFunc"
    return get_sbox(lut, name=name)


