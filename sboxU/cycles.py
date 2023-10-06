#!/usr/bin/sage

from sage.all import Matrix, GF, vector, log, randint, Permutation

from .utils import oplus

from itertools import product as cartesian_prod
from itertools import permutations
from collections import Counter, deque, defaultdict



# !SECTION! Interface between LUTs and `permutations`
# ---------------------------------------------------


def lut_to_permutation(lut):
    """Return the Sage permutation corresponding to the given LUT.
    """
    return Permutation([x + 1 for x in lut])


def permutation_to_lut(perm):
    return [x - 1 for x in perm]



# !SECTION! Cycle decomposition
# -----------------------------


def cycle_decomposition(s):
    """Return a list of all the cycles in the permutation of [0, ...,
    n] with LUT `s`.

    """
    return lut_to_permutation(s).to_cycles()


def cycle_type(s):
    """Return the cycle type of the permutation with LUT `s`, i.e. the
    list of the lengths of all the cycles in the decomposition of `s`.

    """
    result = list(lut_to_permutation(s).cycle_type())
    result.sort(reverse=True)
    return result
    


# !SECTION! Conjugation
# ---------------------

def are_conjugate(a, b):
    """Return True if `a` and `b` are the LUT of two conjugate
    permutations.  Return False otherwise.

    """
    return cycle_type(a) == cycle_type(b)


def nb_conjugates(s):
    """Return the number of conjugates of the permutation with LUT
    `s`.

    """
    sbox_cycle_type = cycle_type(s)
    return SetPartitions(len(s), sbox_cycle_type).cardinality()


def nb_conjugacy_relations(lut):
    """Return the number of conjugacy relations there exists between
    the permutation `s`, and any other permutation conjugated with it.

    """
    sbox_cycle_type = lut_to_perm(lut).cycle_type()
    return product([(cycle_size ** nb_cycles) * factorial(nb_cycles) 
                    for (cycle_size, nb_cycles) in Counter(sbox_cycle_type).items()])



def conjugacy_relations(a, b):
    """An iterator listing all the permutation f such that

    f o a o f^{-1} = b.

    """

    if not are_conjugate(a, b):
        return []

    a_decomp = cycle_decomposition(a)
    b_decomp = cycle_decomposition(b)
    sigma = [0 for x in range(0, len(a))]

    # permutations of the cycles
    for perms in cartesian_prod(*[permutations(range(len(e))) for e in a_decomp]):
        # rotation of the representations of the cycles
        for rotations in cartesian_prod(*[cartesian_prod(range(len(e[0])), repeat=len(e)) for e in a_decomp]):
            # Build the current "conjugating mapping" sigma
            for i_size, cycles_fixed_size in enumerate(a_decomp):
                for i_cycle, cycle in enumerate(cycles_fixed_size):
                    t = deque(b_decomp[i_size][perms[i_size][i_cycle]])
                    t.rotate(rotations[i_size][i_cycle])
                    for i_elt, elt in enumerate(cycle):
                        sigma[elt - 1] = t[i_elt] - 1
            assert comp3(sigma, a, inverse(sigma)) == b
            yield list(sigma)


def conjugates(s):
    """Iterator listing all the conjugates of the permutation `s`."""
    l = len(s)
    sbox_cycle_type = lut_to_permutation(s).cycle_type()
    conjugates = SetPartitions(l, sbox_cycle_type)
    conjugate_v2 = [0 for x in range(l)]
    for conjugate in conjugates:  # conjugate is a Partition
        for c in conjugate:
            c = list(c)
            for i in range(0, len(c)):
                conjugate_v2[c[i] - 1] = c[(i+1) % len(c)] - 1
        # assert are_conjugate(conjugate_v2, lut)
        yield conjugate_v2
