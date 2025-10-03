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
    """The Walsh transform of a function f over of a field of characteristic $c$ is the set

    $W_f(a) = \\sum_x \\rho_c^{ax + f(x)}$

    where ax is a the scalar product of the vectors a and x, where $\\rho_c$ is an order $c$ complex root of unity, and where a takes all possible values. Here, they are ordered using their representation as integers.

    Args:
        s: an S-boxable object

    Returns:
        list: a list `l` such that `l[a]` is equal to $W_f(a)$. Contains integers for vectorial Boolean function, complex numbers otherwise.
    """
    # !CHECK! does the output consists of complex numbers in the case of non binary fields?
    sb = Sb(s)
    if sb.get_output_length() != 1:
        raise Exception("Walsh transform takes as input a boolean function")
    else:
        return cpp_walsh_transform((<S_box>sb).cpp_sb[0])

    
def walsh_spectrum(s):
    """The Walsh spectrum of a function describes the number of time each value appears in its Walsh transform (as returned by `walsh_transform`). For a vectorial function, it counts the number of occurrences of each coefficient in its LAT (as returned by `lat`).

    Args:
        s: an S-boxable object

    Returns:
        Spectrum: a Spectrum instance `d` such that `d[i]` is the number of occurrences of `i` in the LAT of `s`.
    """
    sb = Sb(s)
    result = Spectrum(name="Walsh".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_walsh_spectrum((<S_box>sb).cpp_sb[0],
                           n_threads)
    )
    return result


def absolute_walsh_spectrum(s):
    """The absolute Walsh transform counts the number of occurrences of coefficients with a given absolute value in the Walsh transform for a function (or LAT for a vectorial function).

    Args:
        s: an S-boxable object

    Returns:
        Spectrum: a Spectrum instance `d` such that `d[i]` is the sum of the number of occurrences of `i` and `-i` in the LAT of `s`.
    """
    return walsh_spectrum(s).absolute()
    

def lat(s):
    """The Linear Approximation Table of a function F over a vector space of a field of prime characteristic p is a two dimensional array $W_F$ such that

    $ W_F(a, b) = \\sum_x \\rho_c^{ax + bF(x)} $

    where $ax$ and $bF(x)$ denote the scalar products of $a$ and $x$, and $b$ and $F(x)$; and where $\\rho_c$ is a complex order c root of unity (i.e. $\\rho_c = -1$ for Boolean functions).
    
    Args:
        s: an S-boxable object
    
    Returns:
        list: A list of list `l` such that `l[a][b]` = $W_{F}(a, b)$.
    """
    sb = Sb(s)
    result = cpp_lat((<S_box>sb).cpp_sb[0])
    return result


def invert_lat(l):
    """The LAT is essentially a Fourier transform, meaning that it can be inverted. That is what this function does.

    Args:
        l (list): A list of list corresponding to the LAT of a function.

    Returns:
        An S_box instance whose LAT is `l`.
    """
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
    """
    The XDDT of a function F is a three dimensional array `D` such that `D[a][b] = {x, F(x+a) = F(x)+b}`. The number of rows and columns of the XDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the XDDT of `s`.
    """
    sb=Sb(s)
    result = cpp_xddt((<S_box>sb).cpp_sb[0])
    return result


def yddt(s):
    """
    The YDDT of a function F is a three dimensional array `D` such that `D[a][b] = {F(x), F(x+a) = F(x)+b}`. The number of rows and columns of the YDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the YDDT of `s`.
    """
    sb=Sb(s)
    result = cpp_yddt((<S_box>sb).cpp_sb[0])
    return result


def zddt(s):
    """
    The ZDDT of a function F is a three dimensional array `D` such that `D[a][b] = {x|F(x), F(x+a) = F(x)+b}`. The number of rows and columns of the ZDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the ZDDT of `s`.
    """
    sb=Sb(s)
    result = cpp_zddt((<S_box>sb).cpp_sb[0])
    return result

