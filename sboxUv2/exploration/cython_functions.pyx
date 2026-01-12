# -*- python -*-

from sboxUv2.cython_types cimport *
from sboxUv2.core.sbox import Sb


def non_trivial_sn(mode,s,ne,ns):
    sb = Sb(s)
    result = []
    i = 0
    SW = cpp_non_trivial_sn(<cpp_Integer> mode,(<S_box>sb).cpp_sb[0],<cpp_Integer> ne, <cpp_Integer> ns )
    for sw_u in SW: 
        res_u = []
        for new_s in sw_u:
            new_sb = S_box(name=b"SW-" + sb.name() + b"_" + str(i).encode("UTF-8"))
            new_sb.set_inner_sbox(<cpp_S_box>new_s)
            res_u.append(new_sb)
            i += 1
        result.append(res_u)
    return result

