# -*- python -*-

from sboxUv2.config import N_THREADS


def extract_bases(
        z,
        dimension,
        end_condition="fixed dimension",
        n_threads=N_THREADS):
    """Returns a list containing precisely the GJB basis of each vector space of dimension `d` contained in the list `z`, interpreted like a set of binary vectors.

    The inner working of this algorithms are explained in [AC:BonPerTia19].

    Args:
        - z: a list of positive integers, each being interpreted as a binary vector.
        - dimension: the dimension of the vector spaces to search.
        - end_condition: a string describing what the search is supposed to return. If specified (default is 'fixed dimension', then must be either:
            - 'fixed dimension': the function returns exactly one basis for each space of the given dimension, even when the presence of a space of higher dimension causes multiple lower-dimension spaces to appear.
            - 'just one': the function returns the first vector space of the given dimension found.
            - 'all dimension': if a vector space of a dimension higher than `dimension` is found in `z`, then its basis (which larger than `dimension`) will be in the output.

    Returns:
        A list of list B^i, each list B^i corresponding to the GJB basis of a vector space that is entirely contained in `z`. `z` is always considered to contain 0, even if it is not in it.
    
    """
    return cpp_extract_bases(z, dimension, end_condition, n_threads)


def extract_affine_bases(
        z,
        dimension,
        end_condition="fixed dimension",
        n_threads=N_THREADS):
    """Returns a list containing precisely the GJB basis of each affine space of dimension `d` contained in the list `z`, interpreted like a set of binary vectors. The GJB basis consists of an offset (the smallest element in the affine space), and the GJB basis of its linear part.

    The inner working of this algorithms are explained in [AC:BonPerTia19].

    Args:
        - z: a list of positive integers, each being interpreted as a binary vector.
        - dimension: the dimension of the vector spaces to search.
        - end_condition: a string describing what the search is supposed to return. If specified (default is 'fixed dimension', then must be either:
            - 'fixed dimension': the function returns exactly one basis for each space of the given dimension, even when the presence of a space of higher dimension causes multiple lower-dimension spaces to appear.
            - 'just one': the function returns the first vector space of the given dimension found.
            - 'all dimension': if a vector space of a dimension higher than `dimension` is found in `z`, then its basis (which larger than `dimension`) will be in the output.

    Returns:
        A list of list B^i, each list B^i consisting of an offset followed by the GJB basis of a vector space. The corresponding affine space is entirely contained in `z`.
    
    """
    return cpp_extract_affine_bases(z, dimension, end_condition, n_threads)


