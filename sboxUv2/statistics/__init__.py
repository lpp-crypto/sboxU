# -*- python -*-

"""This module provides tools to study the "statistical" properties of vectorial Boolean functions.

"Statistical" should be understood in the sense of "statistical cryptanalysis", i.e., differential and linear attacks along with all of their variants (boomerang, etc). In practice, this module provides functions that compute the main tables: DDT, LAT, BCT (as well as others), and tools to make sense of their content.

When possible the functions are multi-threaded. In particular, computing "spectra" (i.e., the number of occurrences of the coefficients in a table) is not done by first generating the full table and then counting; instead, rows are generated in parallel, which saves space and increases speed.

"""

from .cython_functions import \
    Spectrum, Spctr, \
    differential_spectrum, ddt, differential_uniformity, is_differential_uniformity_smaller_than, \
    walsh_transform, walsh_spectrum, lat, invert_lat, linearity, \
    boomerang_spectrum, bct, boomerang_uniformity


from .anomalies import \
    ddt_coeff_probability, expected_max_ddt, \
    lat_coeff_probability_permutation, lat_coeff_probability_function, expected_max_lat, expected_max_lat_function, \
    bct_coeff_probability, \
    probability_of_max_and_occurrences, get_proba_func, \
    table_anomaly, table_negative_anomaly
