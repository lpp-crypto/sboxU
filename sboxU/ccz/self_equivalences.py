from sboxU.ccz.affine_equivalence import linear_equivalence
from sboxU.core import F2_trans, is_permutation


#!SECTION! self-equivalence

def self_linear_equivalent_mappings(f):
    """Returns a list of linear permutations L_A,L_B such that L_B o S o L_A = S,
    the result is a list of pairs of F2AffineMap objects.

    """
    return linear_equivalence(f,f,all_mappings=True)

def self_affine_equivalent_mappings(f):
    """Returns the list of the affine permutations A,B such that B o S o A = S.
    """
    if not is_permutation(f):
        raise Exception("argument is not a permutation!")
    result = []
    for cstt_in in f.input_space():
        for cstt_out in f.input_space():
            mappings = linear_equivalence(
                f,
                F2_trans(cstt_out,bit_length=f.get_output_length())*f*F2_trans(cstt_in,bit_length=f.get_input_length()),
                all_mappings=True
            )
            for AB in mappings:
                result.append((cstt_in + AB[0],cstt_out + AB[1]))
    return result

