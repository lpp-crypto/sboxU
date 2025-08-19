
# !SECTION! Threads number

"""The maximum number of threads used in each function. """
MAX_N_THREADS = 8

def n_threads_from_sbox_size(n):
    """Customize this function to adjust the number of threads used by multi-threaded functions.

    The idea is to prevent the use of many threads when each of them actually has very little to do.

    Args:
        - n: the bitlength of an S-box input.

    Returns:
        The number of threads to use. The higher the input, the higher the number of threads.
    """
    if n <= 6:
        return 1
    elif n <= 8:
        return 2
    else:
        return MAX_N_THREADS



# !SECTION! Constants

"""The number of bits of precision to use when performing high precision real arithmetic.

"""
DEFAULT_HIGH_PRECISION = 80



BOLD_START = "\033[1m"
BOLD_END = "\033[0m"
