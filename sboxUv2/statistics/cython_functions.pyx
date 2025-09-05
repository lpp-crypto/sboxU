# -*- python -*-


from sboxUv2.config import n_threads_from_sbox_size

from sboxUv2.core import Sb
from sboxUv2.core cimport *



# !SECTION! Basic tables and spectra

# !SUBSECTION! Differential properties


def differential_spectrum(s):
    sb = Sb(s)
    result = Spectrum(name="Differential".encode("UTF-8"))
    n_threads = n_threads_from_sbox_size(sb.get_input_length())
    result.set_inner_sp(
        cpp_differential_spectrum((<S_box>sb).cpp_sb[0],
                                  n_threads)
    )
    return result


def ddt(s):
    sb = Sb(s)
    result = cpp_ddt((<S_box>sb).cpp_sb[0])
    return result


def differential_uniformity(s):
    sb = Sb(s)
    dif = differential_spectrum(s)
    return dif.maximum()


def is_differential_uniformity_smaller_than(s, u):
    sb = Sb(s)
    return cpp_is_differential_uniformity_smaller_than(
        (<S_box>sb).cpp_sb[0],
        u
    )


# !SUBSECTION! Linear properties


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

        

# !SUBSECTION! Boomerang properties


def boomerang_spectrum(s):
    sb = Sb(s)
    result = Spectrum(name="Boomerang".encode("UTF-8"))
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


# # !SUBSECTION! Degree spectrum 

# def degree_spectrum(s):
#     sb = Sb(s)
#     result = Spectrum(name="Degree".encode("UTF-8"))
#     n_threads = n_threads_from_sbox_size(sb.get_input_length())
#     result.set_inner_sp(
#         cpp_degree_spectrum((<S_box>sb).cpp_sb[0],
#                            n_threads)
#     )
#     return result