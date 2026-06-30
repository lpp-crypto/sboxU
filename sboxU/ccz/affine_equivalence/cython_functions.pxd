# -*- python -*-

from sboxU.cython_types cimport *

from sboxU.core cimport *
from sboxU.statistics cimport *
from sboxU.algorithms cimport *



cdef extern from "../../cpp/ccz/linear_representative.hpp":
    cpp_S_box cpp_le_class_representative(
        const cpp_S_box f,
        cpp_F2AffineMap & A,
        cpp_F2AffineMap & B
    )

    cpp_S_box cpp_le_class_representative(
        const cpp_S_box f,
    )

cdef extern from "../../cpp/ccz/linear_representative.cpp":
    pass


cdef extern from "../../cpp/ccz/linear_equivalence.hpp":
    std_vector[cpp_F2AffineMap] cpp_linear_equivalence_permutations(
        const cpp_S_box f,
        const cpp_S_box g,
        bool all_mappings
    );


cdef extern from "../../cpp/ccz/linear_equivalence.cpp":
    pass
