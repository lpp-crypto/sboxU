# -*- python -*-

from sboxUv2.cython_types cimport *

from sboxUv2.core cimport *
from sboxUv2.statistics cimport *
from sboxUv2.algorithms cimport *



cdef extern from "../../cpp/ccz/linear_representative.hpp":
    cpp_S_box cpp_le_class_representative(
        const cpp_S_box & f,
    )

cdef extern from "../../cpp/ccz/linear_representative.cpp":
    pass


cdef extern from "../../cpp/ccz/linear_equivalence.hpp":
    std_vector[cpp_BinLinearMap] cpp_linear_equivalence_permutations(
        const cpp_S_box f,
        const cpp_S_box g,
        bool all_mappings
    );


cdef extern from "../../cpp/ccz/linear_equivalence.cpp":
    pass
