# -*- python -*-


from sboxUv2.core import Sb, oplus
from sboxUv2.core.sbox import F2_trans
from sboxUv2.config import MAX_N_THREADS
from collections import defaultdict



# !SECTION! XOR equivalence

def xor_equivalence(s, s_prime, all_pairs=True):
    sb, sb_prime = Sb(s), Sb(s_prime)
    if sb.get_input_length() != s_prime.get_input_length():
        return []
    else:
        result = []
        for a in sb.input_space():
            offset = oplus(sb[a], sb_prime[0])
            valid_pair = True
            for x in sb.input_space():
                if (s[oplus(x, a)] != oplus(offset, sb_prime[x])):
                    valid_pair = False
                    break
            if valid_pair:
                result.append([a, offset])
                if not all_pairs:
                    return result
        return result
        

# !SECTION! Linear equivalence

def le_class_representative(s):
    """Computes the smallest member of the linear-equivalence class of the S_boxable object `s`.

    It is assumed that `s` corresponds to a permutation since this function relies on the algorithm presented in [EC:BDCBP03]. A dedicated implementation for the case where s.get_input_length() is at most 8 enables a very computation in this case. Otherwise, defaults to a much slower implementation.

    Args:
        - s: an S_boxable object.

    Returns:
        An `S_box` instance corresponding to the smallest function in the linear equivalence class of `s`, where "smallest" is in the sense of the lexicographic order.
    
    """
    sb = Sb(s)
    if sb.is_invertible():
        result = S_box(name=b"LE(" + sb.name() + b")")
        result.set_inner_sbox(
            # we always try use the "fast" variant of cpp_le_class_representative
            cpp_le_class_representative((<S_box>sb).cpp_sb[0])
        )
        return result
    else:
        raise NotImplemented("Linear representatives can only be computed for permutations")




def linear_equivalence(f, g, all_mappings=False):
    sf = Sb(f)
    sg = Sb(g)
    if len(f) != len(g):
        raise "f and g are of different dimensions!"
    if sf.is_invertible() or sg.is_invertible():
        if sf.is_invertible() and sg.is_invertible():
            return linear_equivalence_permutations(sf, sg, all_mappings=all_mappings)
        else:
            return False # a permutation can only be linear equivalent
                         # to another permutation
    else:
        # !TODO! use Jules table-based algorithm or Itai's algorithm if the degree is maximum
        raise NotImplemented("only permutations are implemented at the moment")
    


def linear_equivalence_permutations(f, g, all_mappings=False):
    """Returns, if it exists, the tuple A, a, B, b where A and B are matrices such that, for all x:

    f = B o g o A

    where "o" denotes functional composition. If no such linear permutations exist, returns an empty list.

    The algorithm used is specified in [EC:BDCBP03].

    """
    sf = Sb(f)
    sg = Sb(g)
    result = cpp_linear_equivalence_permutations(
        (<S_box>sf).cpp_sb[0],
        (<S_box>sg).cpp_sb[0],
        all_mappings
        )
    mappings = []
    for cpp_A in result:
        A = BinLinearMap()
        (<BinLinearMap>A).cpp_blm = new cpp_BinLinearMap()
        (<BinLinearMap>A).cpp_blm[0] = <cpp_BinLinearMap>cpp_A
        mappings.append(A)
    return [(mappings[i], mappings[i+1])
            for i in range(0, len(mappings), 2)]
    
    


# !SECTION! Affine equivalence



def affine_equivalence(f, g):
    sf = Sb(f)
    sg = Sb(g)
    if len(f) != len(g):
        raise "f and g are of different dimensions!"
    if sf.is_invertible() or sg.is_invertible():
        if sf.is_invertible() and sg.is_invertible():
            return affine_equivalence_permutations(sf, sg)
        else:
            return False # a permutation can only be affine equivalent
                         # to another permutation
    else:
        # !TODO! use Jules table-based algorithm in general, and Alain's algorithm for quadratic functions
        raise NotImplemented("only permutations are implemented at the moment")
    

def affine_equivalence_permutations(f, g):
    """Returns, if it exists, the tuple A, a, B, b where A and B are
    matrices and where a and b are integers such that, for all x:

    f(x) = (B o g o A)(x + a) + b,

    where "o" denotes functional composition and "+" denotes XOR. If
    no such affine permutations exist, returns an empty list.

    Internally calls a function written in C++ for speed which returns
    the "Linear Representative" using an algorithm from [EC:BDCBP03].

    """


    sf = Sb(f)
    sg = Sb(g)
    
    if len(f) != len(g):
        raise "f and g are of different dimensions!"
    if not sf.is_invertible():
        raise Exception("first argument is not a permutation!")
    if not sg.is_invertible():
        raise Exception("second argument is not a permutation!")
    table = defaultdict(int)
    n= sf.get_input_length()
    tr = [F2_trans(c,bit_length = n) for c in sf.input_space()] # translations
    for c in sf.input_space():
        table[le_class_representative(tr[c] * sf)] = c
    rs = []
    a = -1
    b = -1    
    for c in sf.input_space():
        g_c = le_class_representative(sg * tr[c])
        if g_c in table.keys():
            a=c
            b = table[g_c]
            rs = g_c
            break
    if a == -1:
        return []
    # !FINISH! 
    l_f = linear_equivalence(tr[b] * sf, rs)
    A_f, B_f = l_f[0][0], l_f[0][1]
    l_g = linear_equivalence(sg * tr[a], rs)
    A_g, B_g = l_g[0][0], l_g[0][1]
    A = A_g.inverse() * A_f
    B = B_f *B_g.inverse()
    a = A.inverse()(a)
    return [A, a, B, b]
