# -*- python -*-

#from sboxUv2.sbox.cython_functions cimport *
from sboxUv2.sbox import Sb
from sboxUv2.config import N_THREADS



def thickness_spectrum(s, spaces=None):
    """
    # !TODO! documentation 
    """
    sb = Sb(s)
    result = Spectrum()
    result.set_inner_sp(
        cpp_thickness_spectrum((<S_box>sb).cpp_sb[0], N_THREADS)
    )
    return result


# !SECTION! Linear equivalence

def le_class_representative(s):
    sb = Sb(s)
    if sb.is_invertible():
        result = S_box(name=b"LE(" + sb.name() + b")")
        result.set_inner_sbox(
            # we always try use the "fast" variant of cpp_le_class_representative
            cpp_le_class_representative((<S_box>sb).cpp_sb[0], 1) 
        )
        return result
    else:
        raise NotImplemented("Linear representatives can only be computed for permutations")
    


# # !SECTION! Affine equivalence



# def affine_equivalence(f, g):
#     sf = Sb(f)
#     sg = Sb(g)
#     if len(f) != len(g):
#         raise "f and g are of different dimensions!"
#     if sf.is_invertible() or sg.is_invertible():
#         if sf.is_invertible() and sg.is_invertible():
#             return affine_equivalence_permutations(sf, sg)
#         else:
#             return False # a permutation can only be affine equivalent
#                          # to another permutation
#     else:
#         raise NotImplemented("only permutations are implemented at the moment")
    

# def affine_equivalence_permutations(f, g):
#     """Returns, if it exists, the tuple A, a, B, b where A and B are
#     matrices and where a and b are integers such that, for all x:

#     f(x) = (B o g o A)(x + a) + b,

#     where "o" denotes functional composition and "+" denotes XOR. If
#     no such affine permutations exist, returns an empty list.

#     Internally calls a function written in C++ for speed which returns
#     the "Linear Representative" using an algorithm from

#     Alex Biryukov, Christophe De Canniere, An Braeken, and Bart
#     Preneel (2003).  "A Toolbox for Cryptanalysis: Linear and Affine
#     Equivalence Algorithms", Advances in Cryptology -- EUROCRYPT 2003,
#     Lecture Notes in Computer Science 2656, E. Biham (ed.),
#     Springer-Verlag, pp. 33--50, 2003.

#     """
#     sf = Sb(f)
#     sg = Sb(g)
#     if len(f) != len(g):
#         raise "f and g are of different dimensions!"
#     if not sf.is_invertible():
#         raise Exception("first argument is not a permutation!")
#     if not sg.is_invertible():
#         raise Exception("second argument is not a permutation!")
#     table = defaultdict(int)
#     tr = [F2_trans(c) for c in sf.input_space()] # translations
#     for b in sf.input_space():
#         table[le_class_representative(tr[c] * f)] = b
#     rs = []
#     a = -1
#     b = -1    
#     for a in sf.input_space():
#         g_c = le_class_representative(g * tr[c])
#         if g_c in table.keys():
#             b = table[g_c]
#             rs = g_c
#             break
#     if a == -1:
#         return []
#     # !CONTINUE!  
#     l_f = linear_equivalence([oplus(f[x], b) for x in range(0, len(f))],
#                              rs)
#     A_f, B_f = l_f[0], l_f[1]
#     l_g = linear_equivalence([g[oplus(x, a)] for x in range(0, len(f))],
#                              rs)
#     A_g, B_g = l_g[0], l_g[1]
#     A = A_g.inverse() * A_f
#     B = B_f * B_g.inverse()
#     a = apply_bin_mat(a, A.inverse())
#     return [A, a, B, b]
