from sage.all import Matrix, log, LinearCode, GF
from sboxU.core import to_bin

# CCZ and EA equivalences ``à la Edel-Pott''
# https://www.yvesedel.de/Papers/nato08.pdf Theorems 9 and 10.


# !TODO! Adapt to handle Sboxes
def ccz_matrix(lut):
    N = int(log(len(lut), 2))
    return Matrix(GF(2), len(lut), 2 * N + 1, 
        [[1] + to_bin((x << N) | lut[x], 2 * N)[::-1] for x in range(0, 2 ** N)]).transpose()


def ea_matrix(lut):
    N = int(log(len(lut), 2))
    return Matrix(GF(2), 2*len(lut)-1, 2 * N + 1,
       [[1] + to_bin((x << N) | lut[x], 2 * N)[::-1] for x in range(0, 2 ** N)] + \
       [to_bin(y, 2 * N + 1)[::-1] for y in range(1, 2 ** N)]).transpose()


def are_code_equivalent(f, g, equivalence="CCZ"):
    if len(f) != len(g):
        raise Exception("f and g are of different sizes!")

    if equivalence == "EA":
        code_f, code_g = LinearCode(ea_matrix(f)), LinearCode(ea_matrix(g))
    elif equivalence == "CCZ":
        code_f, code_g = LinearCode(ccz_matrix(f)), LinearCode(ccz_matrix(g))
    else:
        return None

    res = code_f.is_permutation_equivalent(code_g, algorithm="verbose")

    if res is False:
        return False
    else:
        # !TODO! Handle here the output of the "verbose" code equivalence to obtain the proper affine mapping.
        return True


def are_ea_equivalent_from_code(f, g):
    return are_code_equivalent(f, g, "EA")


def are_ccz_equivalent_from_code(f, g):
    return are_code_equivalent(f, g, "CCZ")
