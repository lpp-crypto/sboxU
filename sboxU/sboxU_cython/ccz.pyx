# -*-python-*- 
# Time-stamp: <2021-08-11 14:47:41 lperrin>

from sboxu_cpp cimport *

def extract_vector(l, a):
	return extract_vector_cpp(l,a)

def extract_bases_fast(l, dimension, word_length, n_threads, end_condition):
	return extract_bases_cpp(l, dimension, word_length, n_threads, end_condition)


def extract_affine_bases_fast(l, dimension, word_length, n_threads, end_condition):
	return extract_affine_bases_cpp(l, dimension, word_length, n_threads, end_condition)	 

def get_lat_zeroes_spaces_fast(l, n, n_threads):
	zeroes = lat_zeroes_cpp(l, n, n_threads)
	return extract_bases_cpp(zeroes, n, 2*n, n_threads, 'fixed dimension'.encode('ascii'))	 


