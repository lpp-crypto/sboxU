# -*- python -*-

"""This module contains the advanced fundamental algorithms that are needed in different parts of sboxU---or can be of independent interest.

More precisely, at this stage, it contains:
- the vector space search ([AC:BonPerTia19]),
- the affine space search easily derived from said vector space search.
- a dedicated data structure that determines the Gauss-Jordan basis of a set of vectors (and thus in particular its rank).

"""

from .cython_functions import \
    extract_bases, extract_affine_bases, \
    BinLinearBasis, is_affine, is_sum_full_rank, \
    complete_basis, complete_basis_reversed, \
    generating_F2AffineMap_r,generating_F2AffineMap, F2AffineMap_from_masks,F2AffineMap_from_range_and_image, \
    F2LinearSystem

