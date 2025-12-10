# -*- python -*-

from sboxUv2.config import MAX_N_THREADS
from sboxUv2.core import oplus
from sboxUv2.core.f2functions import Blm,rank_of_vector_set
from math import log
from random import randint


# !SECTION! Extracting spaces contained in a set

def extract_bases(
        z,
        dimension,
        end_condition="fixed dimension",
        n_threads=MAX_N_THREADS):
    """Returns a list containing precisely the GJB basis of each vector space of dimension `d` contained in the list `z`, interpreted like a set of binary vectors.

    The inner working of this algorithms are explained in [AC:BonPerTia19].

    Args:
        z: a list of positive integers, each being interpreted as a binary vector.
        dimension: the dimension of the vector spaces to search.
        end_condition: a string describing what the search is supposed to return. If specified (default is 'fixed dimension', then must be either:
            - 'fixed dimension': the function returns exactly one basis for each space of the given dimension, even when the presence of a space of higher dimension causes multiple lower-dimension spaces to appear.
            - 'just one': the function returns the first vector space of the given dimension found.
            - 'all dimension': if a vector space of a dimension higher than `dimension` is found in `z`, then its basis (which larger than `dimension`) will be in the output.

    Returns:
        A list of list B^i, each list B^i corresponding to the GJB basis of a vector space that is entirely contained in `z`. `z` is always considered to contain 0, even if it is not in it.
    
    """
    return cpp_extract_bases(z,
                             dimension,
                             n_threads,
                             end_condition.encode("UTF-8"))
                             


def extract_affine_bases(
        z,
        dimension,
        end_condition="fixed dimension",
        n_threads=MAX_N_THREADS):
    """Returns a list containing precisely the GJB basis of each affine space of dimension `d` contained in the list `z`, interpreted like a set of binary vectors. The GJB basis consists of an offset (the smallest element in the affine space), and the GJB basis of its linear part.

    The inner working of this algorithms are explained in [AC:BonPerTia19].

    Args:
        z: a list of positive integers, each being interpreted as a binary vector.
        dimension: the dimension of the vector spaces to search.
        end_condition: a string describing what the search is supposed to return. If specified (default is 'fixed dimension', then must be either:
            - 'fixed dimension': the function returns exactly one basis for each space of the given dimension, even when the presence of a space of higher dimension causes multiple lower-dimension spaces to appear.
            - 'just one': the function returns the first vector space of the given dimension found.
            - 'all dimension': if a vector space of a dimension higher than `dimension` is found in `z`, then its basis (which larger than `dimension`) will be in the output.

    Returns:
        A list of list B^i, each list B^i consisting of an offset followed by the GJB basis of a vector space. The corresponding affine space is entirely contained in `z`.
    
    """
    return cpp_extract_affine_bases(z,
                                    dimension,
                                    n_threads,
                                    end_condition.encode("UTF-8"))



# !SECTION! The BinLinearBasis class


# !SUBSECTION! The class itself 

cdef class BinLinearBasis:
    # !TODO! documentation of the BinLinearBasis class 
    def __init__(self, std_vector[BinWord] l):
        self.cpp_lb = new cpp_BinLinearBasis(l)


    # !TODO! destructor for BinLinearBasis 
    # def __dealloc__(self):
    #     self.cpp_lb[0].destruct()
    #     free(self.cpp_lb)

        
    def __iter__(self) -> BinWord:
        for x in self.cpp_lb[0].get_basis():
            yield x

            
    def __len__(self) -> int:
        return self.cpp_lb[0].rank()


    def __str__(self) -> str:
        result = "("
        for x in self.cpp_lb[0].get_basis():
            result += "{:x}, ".format(x)
        return result[:-2] + ")"

    
    def basis_vectors(self) -> list:
        return [x for x in self.cpp_lb[0].get_basis()]

    
    def rank(self) -> int:
        return self.cpp_lb[0].rank()

    
    def add_to_span(self, BinWord x) -> bool:
        return self.cpp_lb[0].add_to_span(x)
        

    def span(self) -> std_vector[BinWord]:
        return self.cpp_lb[0].span()
    

    def is_in_span(self, BinWord x) -> bool:
        return self.cpp_lb[0].is_in_span(x)

    def __eq__(self, BinLinearBasis b) -> bool:
        return self.basis_vectors()==b.basis_vectors()



# !SUBSECTION! Using BinLinearBasis 

def is_affine(l,give_basis=False):
    b = BinLinearBasis([])
    V = [oplus(l[0], x) for x in l] 
    n = round(log(len(V), 2))
    if 2**n != len(V):
        if give_basis:
            return False,None,None
        else :
            return False
    for x in V:
        b.add_to_span(x)
    if give_basis : 
        if b.rank() != n:
            return False,None,None
        else:
            return True,l[0],[v for v in b]
    else :
        return b.rank()==n


def complete_basis(basis,N):
    if rank_of_vector_set(basis)!=len(basis): 
        raise Exception("in complete_basis: the input must be independent! input={}".format(basis))
    r = len(basis)
    e_i = 1
    res=basis
    while r < N:
        new_basis = res + [e_i]
        new_r = rank_of_vector_set(new_basis)
        if new_r > r:
            res = new_basis[:]
            r = new_r
        e_i += 1
    return res

def complete_basis_reversed(basis,N): 
    if rank_of_vector_set(basis) != len(basis):
        raise Exception("in complete_basis: the input must be independent! input={}".format(basis))
    r = len(basis)
    e_i = 1
    while r < N:
        new_basis = [e_i] + basis
        new_r = rank_of_vector_set(new_basis)
        if new_r > r:
            basis = new_basis[:]
            r = new_r
        e_i += 1
    return basis

def generating_BinLinearMap_r(basis,N)-> BinLinearMap:
    return Blm(complete_basis_reversed(basis,N),N,N)

def generating_BinLinearMap(basis, N)-> BinLinearMap:
    """Returns a BinLinearMap L from N bits to N bits such that L(1 << i) = basis[i] for
    all i < len(basis) and such that L is a bijection .
    """
    return Blm(complete_basis(basis,N),N,N)

def BinLinearMap_from_masks(masks, N,M)-> BinLinearMap:
    """Returns a BinLinearMap L from N bits to M bits such that L(1 << i) = basis[i] for
    all i < len(basis) and L(1<<i)=0 for all i >= len(basis).
    """
    return Blm(masks + [0 for _ in range(len(masks),N)],N,M)

def BinLinearMap_from_range_and_image(inputs,outputs,N,M)-> BinLinearMap:
    return BinLinearMap_from_masks(outputs,M,N)*(generating_BinLinearMap(inputs,M).inverse())
    

# !SECTION!  Solving Linear Systems


cdef class F2LinearSystem:
    """Implements a linear system of equation over F_2.

    Is optimized for the specific case of large systems where equations are added one by one. The rank is evaluated in real time (i.e., every time an equation is added), while the solutions are only obtained once the system is fully known.

    """

    def __init__(self, n_variables: int):
        self.cpp_ls = new cpp_F2LinearSystem(n_variables)


    # !TODO! destructor for F2LinearSystem 
    # def __dealloc__(self):
    #     pass

    def add_equation(self, variable_indices):
        self.cpp_ls[0].add_equation(<std_vector[unsigned int]>variable_indices)

    def remove_solution(self, sol):
        self.cpp_ls[0].remove_solution(<std_vector[unsigned int]>sol)

    def rank(self):
        return self.cpp_ls[0].rank()

    def kernel(self):
        return self.cpp_ls[0].kernel_as_BinWords()
    
    def __str__(self) -> str:
        return self.cpp_ls[0].to_string().decode("UTF-8")
