# -*- python -*-

"""This module contains the advanced fundamental algorithms that are needed in different parts of sboxU---or can be of independent interest.

More precisely, at this stage, it contains:
- the vector space search ([AC:BonPerTia19]),
- the affine space search easily derived from said vector space search.
- a dedicated data structure that determines the Gauss-Jordan basis of a set of vectors (and thus in particular its rank).

"""

from .cython_functions import \
    extract_bases, extract_affine_bases, \
    BinLinearBasis, is_affine, \
    complete_basis, complete_basis_reversed, \
    generating_BinLinearMap_r,generating_BinLinearMap, BinLinearMap_from_masks,BinLinearMap_from_range_and_image, \
    F2LinearSystem,BinLinearBigBasis


from .invariants import all_invariants_up_to_degree