# -*-python-*- 
# Time-stamp: <2021-08-12 17:25:04 lperrin>

from sboxu_cpp cimport *
from math import log


BIG_SBOX_THRESHOLD = 128
DEFAULT_N_THREADS  = 2


def extract_vector(l, a):
    return extract_vector_cpp(l,a)


def extract_bases_fast(l, dimension, word_length, n_threads, end_condition):
    return extract_bases_cpp(l, dimension, word_length, n_threads, end_condition)


def extract_affine_bases_fast(l, dimension, word_length, n_threads, end_condition):
    return extract_affine_bases_cpp(l, dimension, word_length, n_threads, end_condition)


def get_lat_zeroes_spaces_fast(l, n, n_threads):
    zeroes = lat_zeroes_cpp(l, n, n_threads)
    return extract_bases_cpp(zeroes, n, 2*n, n_threads, 'fixed dimension'.encode('ascii'))


def sigma_multiplicities(s, k, n_threads=None):
    """The multiset \\Sigma_F^k(0) as defined in
    https://seta-2020.org/assets/files/program/papers/paper-44.pdf

    """
    if n_threads == None:
        if len(s) > BIG_SBOX_THRESHOLD:
            n_threads = DEFAULT_N_THREADS
        else:
            n_threads = 1
    n = int(log(len(s), 2))
    return sigma_multiplicities_cpp(s, k, n, n_threads)
