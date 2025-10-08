# -*- python -*-

from libcpp.vector cimport vector as std_vector
from libcpp.string cimport string
from libcpp cimport bool
from libc.stdint cimport uint64_t as BinWord
from libc.stdint cimport int64_t, uint8_t
from libcpp.utility cimport pair

ctypedef std_vector[int64_t] FpWord
ctypedef std_vector[uint8_t] Bytearray
ctypedef int64_t cpp_Integer
