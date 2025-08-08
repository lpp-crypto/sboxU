# -*- python -*-

"""This module provides tools to study the "statistical" properties of vectorial Boolean functions.

"Statistical" should be understood in the sense of "statistical cryptanalysis", i.e., differential and linear attacks along with all of their variants (boomerang, etc). In practice, this module provides functions that compute the main tables: DDT, LAT, BCT (as well as others), and tools to make sense of their content.

When possible the functions are multi-threaded. In particular, computing "spectra" (i.e., the number of occurrences of the coefficients in a table) is not done by first generating the full table and then counting; instead, rows are generated in parallel, which saves space and increases speed.

"""

from .cython_functions import *

from .anomalies import *
