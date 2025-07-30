# Created by Jules Baudrin on 06/07/2025
# Modified by Jules Baudrin on 30/07/2025

import os
from .sboxu_cpp cimport *

def is_linearly_self_equivalent_from_lat(lat, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """ Checks whether a function F is linearly self-equivalent (*) from its linear approximation table (LAT).
     (*) i.e. if it exists a non-trivial pair of linear bijections (A, B) satisfying B o F o A = F.

    More information on the search algorithms can be found in the corresponding C++ code.

    Args:
        lat: The linear approximation table of F.
        algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                    - "alt_partition_diag_mappings" (default)
                    - "alt_partition"
                    - "std_partition_diag_mappings"
        number_of_threads (int): The number of threads to use. By default, it is set to 8.

    Returns:
        A pair of look-up tables corresponding to a non-trivial (A, B) satisfying B o F o A = F, *if it exists*.
        Return False otherwise.

    """
    res = cpp_is_linearly_self_equivalent_from_lat(lat, algorithm.encode('ascii'), number_of_threads)
    if res == ([], []):
        return False
    else:
        return res

def linear_automorphisms_from_lat(lat, algorithm="alt_partition_diag_mappings", number_of_threads=8):
    """ Compute the linear automorphisms (*) of a function F from its linear approximation table (LAT).
        (*) i.e. all pairs of linear bijections (A, B) satisfying B o F o A = F.

        More information on the search algorithms can be found in the corresponding C++ code.

        Args:
            lat: The linear approximation table of F.
            algorithm (str): A string describing the search algorithm to use. algorithm must be chosen among:
                        - "alt_partition_diag_mappings" (default)
                        - "alt_partition"
                        - "std_partition_diag_mappings"
            number_of_threads (int): The number of threads to use. By default, it is set to 8.

        Returns:
            list: A list of pairs of look-up tables corresponding to all pairs of linear bijections (A, B) satisfying B o F o A = F.

    """
    return cpp_linear_automorphisms_from_lat(lat, algorithm.encode('ascii'), number_of_threads)