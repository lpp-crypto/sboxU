# -*- python -*-

"""This module contains the advanced fundamental algorithms that are needed in different parts of sboxU---or can be of independent interest.

More precisely, at this stage, it contains:
- the vector space search ([AC:BonPerTia19]),
- the affine space search easily derived from said vector space search.
- a dedicated data structure that determines the Gauss-Jordan basis of a set of vectors (and thus in particular its rank).

"""

from .cython_functions import \
    extract_bases, extract_affine_bases, \
    Linear_basis

