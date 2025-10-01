# -*- python -*-


from sboxUv2.config import n_threads_from_sbox_size

from sboxUv2.core import Sb
from sboxUv2.core cimport *




# !SECTION! Differential properties


def differential_spectrum(s):
    """The differential spectrum of an S-box counts the number of entries in the DDT that are equal to each value.

    Args:
        s: An S-boxable object.
    
    Returns:
        Spectrum: A `Spectrum` instance `d` such that `d[k]` is equal to the number of occurrences of the coefficient `k` in the DDT of `s`.
    
    """
    sb = Sb(s)
    result = Spectrum(name="Differential".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_differential_spectrum((<S_box>sb).cpp_sb[0],
                                  n_threads)
    )
    return result


def ddt(s):
    """The Difference Distribution Table of a function F is a two dimensional array `D` such that `D[a][b] = #{x, F(x+a) = F(x)+b}`. The number of rows and columns of the DDT depends on the dimensions of the input and output of the S-box (respectively). Depending on the underlying arithmetic of `s`, `+` corresponds either to a XOR (in F_2) or to a modular addition.

    Args:
        s: An S-boxable object.

    Returns:
        list: A list of lists corresponding to a 2-dimensional array containing the DDT of `s`.
    
    """
    sb = Sb(s)
    result = cpp_ddt((<S_box>sb).cpp_sb[0])
    return result


def differential_uniformity(s):
    """The differential uniformity of a function F is the maximum over all a != 0 and all b of the number of solutions x of the equation F(x+a)=F(x)+b.

    Args:
        s: An S-boxable object.

    Returns:
        int: The differential uniformity of the function corresponding to `s`.
    """
    sb = Sb(s)
    dif = differential_spectrum(s)
    return dif.maximum()


def is_differential_uniformity_smaller_than(s, u):
    """Tests whether the differential uniformity of a function is below or equal to a given threshold. If the answer is "no", then its execution can be much smaller than a proper computation of the differential uniformity.

    Args:
        s: An S-boxable object.
        u (int): The threshold such differential uniformity <= u.

    Returns:
        bool: True if and only if the differential uniformity of the function corresponding to `s` is at most equal to `u`.
    
    """
    sb = Sb(s)
    return cpp_is_differential_uniformity_smaller_than(
        (<S_box>sb).cpp_sb[0],
        u
    )


# !SECTION! Linear properties


def walsh_transform(s):
    sb = Sb(s)
    if sb.get_output_length() != 1:
        raise Exception("Walsh transform takes as input a boolean function")
    else:
        return cpp_walsh_transform((<S_box>sb).cpp_sb[0])

    
def walsh_spectrum(s):
    sb = Sb(s)
    result = Spectrum(name="Walsh".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_walsh_spectrum((<S_box>sb).cpp_sb[0],
                           n_threads)
    )
    return result


def absolute_walsh_spectrum(s):
    return walsh_spectrum(s).absolute()
    

def lat(s):
    sb = Sb(s)
    result = cpp_lat((<S_box>sb).cpp_sb[0])
    return result


def invert_lat(l):
    sb = S_box(name="LAT^-1")
    (<S_box>sb).set_inner_sbox(<cpp_S_box>cpp_invert_lat(l))
    return sb


def linearity(s):
    sb = Sb(s)
    wal = walsh_spectrum(s)
    return wal.maximum()

        

# !SECTION! Boomerang properties

# !SUBSECTION! BCT 

def boomerang_spectrum(s):
    sb = Sb(s)
    result = Spectrum(name="BCT".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_boomerang_spectrum((<S_box>sb).cpp_sb[0],
                               n_threads)
    )
    return result


def bct(s):
    sb = Sb(s)
    result = cpp_bct((<S_box>sb).cpp_sb[0])
    return result


def boomerang_uniformity(s):
    sb = Sb(s)
    bom = boomerang_spectrum(s)
    return bom.maximum()



# !SUBSECTION! FBCT

def fbct_spectrum(s):
    sb = Sb(s)
    result = Spectrum(name="F-BCT".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_fbct_spectrum((<S_box>sb).cpp_sb[0],
                          n_threads)
    )
    return result



def fbct(s):
    sb = Sb(s)
    result = cpp_fbct((<S_box>sb).cpp_sb[0])
    return result

# SECTION xddt and co

def xddt(s):
    sb=Sb(s)
    result = cpp_xddt((<S_box>sb).cpp_sb[0])
    return result

