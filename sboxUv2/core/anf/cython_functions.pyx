# cython_functions.pyx


# # Wrapper Python pour degree component
# def degree_component(cpp_S_box f):
#     return cpp_degree_component(f)

# # Wrapper Python pour monomial degree spectrum
# def monomial_degree_spectrum_component(cpp_S_box f):
#     return cpp_monomial_degree_spectrum_component(f)

# Import des wrappers Python de core/sbox
from ..sbox.cython_functions cimport S_box, cpp_S_box

# Import des wrappers Python de core/statistics (Spectrum)
from ...statistics.cython_functions cimport Spectrum, cpp_Spectrum


from .cython_functions cimport cpp_degree_spectrum

def degree_spectrum(S_box s):
    py_result = Spectrum(name=b"Degree")  
    py_result.set_inner_sp(
        cpp_degree_spectrum(<const cpp_S_box&>s.cpp_sb[0])
    )
    return py_result
