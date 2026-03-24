# -*- python -*-


from sboxU.config import n_threads_from_sbox_size

from sboxU.core import get_sbox
from sboxU.core cimport *


from cython.operator cimport dereference



# !SECTION! Differential properties


def differential_spectrum(s):
    """The differential spectrum of an S-box counts the number of entries in the DDT that are equal to each value.
    
    This function does not store the DDT in memory before counting, and uses openMP multi-threading to speed things further.

    Args:
        s: An S-boxable object.
    
    Returns:
        Spectrum: A `Spectrum` instance `d` such that `d[k]` is equal to the number of occurrences of the coefficient `k` in the DDT of `s`.
    
    """
    sb = get_sbox(s)
    result = Spectrum(name="Differential".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_differential_spectrum(dereference((<S_box>sb).cpp_sb),
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
    sb = get_sbox(s)
    result = cpp_ddt(dereference((<S_box>sb).cpp_sb))
    return result


def differential_uniformity(s):
    """The differential uniformity of a function F is the maximum over all a != 0 and all b of the number of solutions x of the equation F(x+a)=F(x)+b.

    Args:
        s: An S-boxable object.

    Returns:
        int: The differential uniformity of the function corresponding to `s`.
    """
    sb = get_sbox(s)
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
    sb = get_sbox(s)
    return cpp_is_differential_uniformity_smaller_than(
        dereference((<S_box>sb).cpp_sb),
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
    sb = get_sbox(s)
    if sb.get_output_length() != 1:
        raise Exception("Walsh transform takes as input a boolean function")
    else:
        return cpp_walsh_transform(dereference((<S_box>sb).cpp_sb))

    
def walsh_spectrum(s):
    """The Walsh spectrum of a function describes the number of time each value appears in its Walsh transform (as returned by `walsh_transform`). For a vectorial function, it counts the number of occurrences of each coefficient in its LAT (as returned by `lat`).
    
    This function does not store the LAT in memory before counting, and uses openMP multi-threading to speed things further.

    Args:
        s: an S-boxable object

    Returns:
        Spectrum: a Spectrum instance `d` such that `d[i]` is the number of occurrences of `i` in the LAT of `s`.
    """
    sb = get_sbox(s)
    result = Spectrum(name="Walsh".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_walsh_spectrum(dereference((<S_box>sb).cpp_sb),
                           n_threads)
    )
    return result


def absolute_walsh_spectrum(s):
    """The absolute Walsh transform counts the number of occurrences of coefficients with a given absolute value in the Walsh transform for a function (or LAT for a vectorial function).
    
    This function does not store the LAT in memory before counting, and uses openMP multi-threading to speed things further.

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
    sb = get_sbox(s)
    result = cpp_lat(dereference((<S_box>sb).cpp_sb))
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
    """The linearity is the maximum module/absolute value to be found among non-trivial coefficients in the Walsh transform/LAT of a function.

    Args:
        s: an S-boxable object

    Returns:
        number: In F_2, an integer; in F_p, a real number.
    """
    sb = get_sbox(s)
    wal = walsh_spectrum(s)
    return wal.maximum()

        

# !SECTION! Boomerang properties

# !SUBSECTION! BCT 

def boomerang_spectrum(s):
    """The boomerang spectrum captures the number of occurrences of each coefficient in the BCT of a vectorial Boolean function (as returned by the `bct` function).

    This function does not store the full BCT in memory before counting, and uses openMP multi-threading to speed things further.

    Args:
        s: an S-boxable object

    Returns:
        Spectrum: a Spectrum instance `d` such that `d[i]` is the number of occurrences of `i` in the BCT of `s`.
    """
    sb = get_sbox(s)
    result = Spectrum(name="BCT".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_boomerang_spectrum(dereference((<S_box>sb).cpp_sb),
                               n_threads)
    )
    return result


def bct(s):
    """The Boomerang Connectivity Table was introduced in [EC:CHPS+18] to better capture what happens at the transition between the forward and the backward trail in a boomerang attack against an SPN.

    It is a two dimensional array B such that, for all a and b:

    $B(a,b) = \\#\\{x, S^{-1}(S(x)+b) + S^{-1}(S(x+a)+b)=a\\}$
    """
    sb = get_sbox(s)
    result = cpp_bct(dereference((<S_box>sb).cpp_sb))
    return result


def boomerang_uniformity(s):
    sb = get_sbox(s)
    bom = boomerang_spectrum(s)
    return bom.maximum()



# !SUBSECTION! FBCT

def fbct_spectrum(s):
    sb = get_sbox(s)
    result = Spectrum(name="F-BCT".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_fbct_spectrum(dereference((<S_box>sb).cpp_sb),
                          n_threads)
    )
    return result



def fbct(s):
    sb = get_sbox(s)
    result = cpp_fbct(dereference((<S_box>sb).cpp_sb))
    return result


# !SECTION! Differential-Linear properties


def dlct(s):
    """The Differential-Linear Connectivity Table (DLCT) of a vectorial Boolean function F,
<<<<<<< HEAD
    as introduced in [AC:BCCD+19] (Cid, Huang, Peyrin, Sasaki, Song).
=======
    as introduced in [EC:BDKW19] (Bar-On, Dunkelman, Keller, Weizman).
>>>>>>> 7a033b0 (feat: implement Differential-Linear Connectivity Table (DLCT))

    The DLCT is defined as the 2^n x 2^n table:

        DLCT[a][b] = sum_{x in F_2^n} (-1)^{b · (F(x) XOR F(x XOR a))}

    Equivalently, each row of the DLCT is the Fast Walsh-Hadamard Transform of the
    corresponding row of the DDT.  For a != 0 and b != 0, DLCT[a][b] / 2^n measures
    the correlation of the differential transition a -> * with the output linear mask b,
    which governs the success probability of differential-linear distinguishers.

    Args:
        s: an S-boxable object over F_2.

    Returns:
        list: A list of 2^n lists of 2^n integers giving the full DLCT.

    """
    sb = get_sbox(s)
    n = sb.get_input_length()
    N = 1 << n
    lut = [sb(x) for x in range(N)]
    table = []
    for a in range(N):
        # Build DDT row: row[v] = #{x : F(x) XOR F(x XOR a) = v}
        row = [0] * N
        for x in range(N):
            row[lut[x] ^ lut[x ^ a]] += 1
        # Apply in-place Fast Walsh-Hadamard Transform to get DLCT row
        h = 1
        while h < N:
            for i in range(0, N, h << 1):
                for j in range(i, i + h):
                    u, v = row[j], row[j + h]
                    row[j] = u + v
                    row[j + h] = u - v
            h <<= 1
        table.append(row)
    return table


def dlct_spectrum(s):
    """The DLCT spectrum counts the occurrences of each value in the Differential-Linear
    Connectivity Table for non-trivial pairs (a != 0, b != 0).

    Args:
        s: an S-boxable object over F_2.

    Returns:
        Spectrum: a Spectrum instance `d` such that `d[v]` is the number of pairs (a, b) with
        a != 0 and b != 0 for which DLCT[a][b] = v.

    """
    table = dlct(s)
    N = len(table)
    result = Spectrum(name="DLCT".encode("UTF-8"))
    result.incr_by_counting([table[a][b] for a in range(1, N) for b in range(1, N)])
    return result


def dlct_uniformity(s):
    """The DLCT uniformity is the maximum absolute value of DLCT[a][b] over all a != 0
    and b != 0.  A lower value indicates stronger resistance to differential-linear attacks.

    Args:
        s: an S-boxable object over F_2.

    Returns:
        int: the DLCT uniformity of the function corresponding to `s`.

    """
    table = dlct(s)
    N = len(table)
    # Cannot use dlct_spectrum(s).maximum() here: Spectrum.maximum() only
    # considers positive keys, so it would miss cases where the largest
    # absolute value is negative.
    return max(abs(table[a][b]) for a in range(1, N) for b in range(1, N))


# SECTION xddt and co

def xddt(s):
    """
    The XDDT of a function F is a three dimensional array `D` such that `D[a][b] = {x, F(x+a) = F(x)+b}`. The number of rows and columns of the XDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the XDDT of `s`.
    """
    sb=get_sbox(s)
    result = cpp_xddt(dereference((<S_box>sb).cpp_sb))
    return result


def yddt(s):
    """
    The YDDT of a function F is a three dimensional array `D` such that `D[a][b] = {F(x), F(x+a) = F(x)+b}`. The number of rows and columns of the YDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the YDDT of `s`.
    """
    sb=get_sbox(s)
    result = cpp_yddt(dereference((<S_box>sb).cpp_sb))
    return result


def zddt(s):
    """
    The ZDDT of a function F is a three dimensional array `D` such that `D[a][b] = {x|F(x), F(x+a) = F(x)+b}`. The number of rows and columns of the ZDDT depends on the dimensions of the input and output of the S-box (respectively).

    Args :
        s: an S-boxable object over F_2
    
    Returns :
        list: A list of lists of lists corresponding to a 3-dimensional array containing the ZDDT of `s`.
    """
    sb=get_sbox(s)
    result = cpp_zddt(dereference((<S_box>sb).cpp_sb))
    return result

# !SECTION! Linear structures

def linear_structures(s):
    """
    Args :
        s: an S-boxable object over F_2 of output length 1
    
    Returns: 
        list,list : A pair of lists `[l_0, l_1]` such that, for all `a` in l_e (with e in [0,1]),
        f(x+a) + f(x) = e,
        for all x (where `+` corresponds to a XOR). 0 does not appear in `l_0` as it would always be present.
    """
    sb=get_sbox(s)
    result=cpp_linear_structures(dereference((<S_box>sb).cpp_sb))
    return result

def linear_structures_vectorial(s):
    """
    Args :
        s: an S-boxable object over F_2 
    
    Returns : 
        dict : A dictionnary `d`, where `d[c]` is a pair of lists `[l_0,l_1]` such that, for all `a` in l_e (with e in [0,1]),
        c. (s(x+a) + s(x)) = e,
        for all x (where `+` corresponds to a XOR), where `.` is the
        scalar product.
        The dictionnary only contains keys where `l_0` or `l_1` are
        non-trivial, i.e. where `l_0` contains more than just 0 and/or
        where `l_1` is non-empty.

    """
    sb=get_sbox(s)
    result=cpp_linear_structures_vectorial(dereference((<S_box>sb).cpp_sb))
    return result

def linear_structures_vectorial_spectrum(s):
    """
    Args :
        s: an S-boxable object over F_2 
    
    Returns : 
        Spectrum: a Spectrum instance `d` such that `d[c]` is the number of linear structures of the Boolean function c.s , where `.` is the
        scalar product. The spectrum only contains keys for which the number of linear structures is non-zero.

    """
    sb=get_sbox(s)
    result=Spectrum(name="Linear Structures".encode("UTF-8"))
    result.set_inner_sp(cpp_linear_structures_vectorial_spectrum(dereference((<S_box>sb).cpp_sb)))
    return result
