# -*- python -*-

from sboxUv2.config import MAX_N_THREADS
from sboxUv2.core import oplus
from sboxUv2.core.f2functions import Blm,rank_of_vector_set
from math import log
from random import randint

from cython.operator cimport dereference



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
    """Stores the Gauss-Jordan basis of a vector space.

    For a given vector space, this basis is uniquely defined (it is for example the smallest one if you were to rank all the basis of a vector space in lexicographic order), meaning that testing the equality of the BinLinearBasis of two vector spaces is enough to test their equality.

    It provides efficient method to add a new dimension to the span of the corresponding vector space, and to test if an element is the corresponding vector space.

    It is a thin wrapper for the `cpp_BinLinearBasis` C++ class. Check out its documentation for more information.

    It only supports vectors of length at most 64 since, under the hood, it stores each vector in a `uint64_t`.
    
    """
    
    def __init__(self, std_vector[BinWord] l):
        """Creates a BinLinearBasis instance that corresponds to the smallest vector space containing all the vectors in `l`.

        Args:
            l: a list of integers interpreted as a list of binary vectors.
        
        """
        self.cpp_lb = make_unique[cpp_BinLinearBasis](<std_vector[BinWord]> l)

        
    def __iter__(self) -> BinWord:
        """An iterator for the basis elements.

        For example, for the following snippet::

        >>> print([x  for x in BinLinearBasis([1, 2])])
        [1, 2]

         3 will not appear.

        Yields:
            The next element in the basis.
        
        """
        for x in dereference(self.cpp_lb).get_basis():
            yield x

            
    def __len__(self) -> int:
        """Returns:
            The number of elements in the basis, i.e., the rank of the set of vectors constituting this BinLinearBasis.
        
        """
        return dereference(self.cpp_lb).rank()


    def __add__(self, b : BinLinearBasis | list[BinWord]) -> BinLinearBasis:
        """Combine two BinLinearBasis instances.

        Args:
            b: either a BinLinearBasis, or a list of BinWord that corresponds to the content of a BinLinearBasis.

        Returns:
            A BinLinearBasis spanning the content of both of the inputs bases. The dimension of the result is at least the maximum of the dimensions of the inputs, and at most the sum of their dimensions.
        """
        if isinstance(b, BinLinearBasis):
            x = b
        elif isinstance(b, list):
            x = BinLinearBasis(b)
        else:
            raise NotImplemented
        return BinLinearBasis(
            dereference((<BinLinearBasis>self).cpp_lb).add(dereference((<BinLinearBasis>x).cpp_lb)).get_basis()
        )
        
    
    def rank(self) -> int:
        """Returns:
            The rank of the set of vectors constituting this BinLinearBasis, i.e., the number of elements in the basis.
        
        """
        return dereference(self.cpp_lb).rank()


    def basis_vectors(self) -> list[BinWord]:
        """Returns:
            An ordered list containing the integers corresponding to the vectors forming the basis of the vector space corresponding to this BinLinearBasis instance.
        
        """
        return [x for x in dereference(self.cpp_lb).get_basis()]

    
    def add_to_span(self, BinWord x) -> bool:
        """Adds a vector to the vector space spanned by this BinLinearBasis. If `x` is actually already in said space, i.e., if it is spanned by the current basis vectors, then does nothing. Otherwise, builds a new echelonized basis that will contain one more element, this new basis spanning the previous space V as well as x+V.

        The echelonizing step means that this function has a run time that is proportional to the current size of the basis.

        Args:
            x: The integer representing the binary vector to add to the vector space corresponding to this BinLinearBasis.
        
        Returns:
            True if adding the element has indeed increased the dimension of the vector space, False if the BinLinearBasis is actually left unchanged by this operation (i.e. if the element was already in the vector space it spans. This output is equal to `self.is_in_span(x)`.
        
        """
        return dereference(self.cpp_lb).add_to_span(x)
        

    def span(self) -> std_vector[BinWord]:
        """Computes the full vector space spanned by the vectors in the basis.

        Returns:
            An ordered list of integers, one per vector in the vector space spanned by this basis.
        """
        return dereference(self.cpp_lb).span()
    

    def is_in_span(self, BinWord x) -> bool:
        """Tests if a vector is contained in the vector space spanned by this basis.

        Args:
            x: the integer representation of the binary vector to test.

        Returns:
            True if the binary vector corresponding to `x` is spanned by the vectors in this basis (or, equivalently, if it is in the vector space that has this basis as its Guss-Jordan basis), False otherwise.
        """
        return dereference(self.cpp_lb).is_in_span(x)

    
    def __eq__(self, BinLinearBasis b) -> bool:
        """An equality test between two BinLinearBasis.

        Args:
            b: a BinLinearBasis to test

        Returns:
            True if `b` contains the same vectors as this BinLinearBasis instance, False otherwise. Equivalently, returns True if and only if the vector spaces spanned by `b` and by this basis are the same.
        """
        return self.basis_vectors()==b.basis_vectors()


    def __str__(self) -> str:
        result = "("
        for x in dereference(self.cpp_lb).get_basis():
            result += "{:x}, ".format(x)
        return result[:-2] + ")"

    

# !SUBSECTION! Using BinLinearBasis

def is_sum_full_rank(b0: list[BinWord] | BinLinearBasis, b1: list[BinWord] | BinLinearBasis) -> bool:
    # checking b0
    if isinstance(b0, list):
        x0 = BinLinearBasis(<std_vector[BinWord]> b0)
    elif isinstance(b0, BinLinearBasis):
        x0 = b0
    else:
        raise Exception("wrong input type for is_sum_full_rank, expected list[BinWord] or BinLinearBasis, got {}".format(type(b0)))
    # checking b1
    if isinstance(b1, list):
        x1 = BinLinearBasis(<std_vector[BinWord]> b1)
    elif isinstance(b1, BinLinearBasis):
        x1 = b1
    else:
        raise Exception("wrong input type for is_sum_full_rank, expected list[BinWord] or BinLinearBasis, got {}".format(type(b1)))
    # actual call
    return cpp_is_sum_full_rank(
        dereference((<BinLinearBasis>x0).cpp_lb),
        dereference((<BinLinearBasis>x1).cpp_lb)
    )
    

def is_affine(l: list[BinWord], give_basis=False):
    # !TODO! docstring for algorithm.is_affine
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


def complete_basis(basis: list[BinWord], N: int) -> list[BinWord]:
    # !TODO! docstring for algorithm.complete_basis
    # !TODO! replace cython complete_basis by a wrapper for cpp_complete_basis 
    if rank_of_vector_set(basis)!=len(basis): 
        raise Exception("in complete_basis: the inputs must be independent! input={}".format(basis))
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


def complete_basis_reversed(basis: list[BinWord], N: int) -> list[BinWord]: 
    # !TODO! docstring for algorithm.complete_basis_reversed
    if rank_of_vector_set(basis) != len(basis):
        raise Exception("in complete_basis: the inputs must be independent! input={}".format(basis))
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


def generating_BinLinearMap_r(basis: list[BinWord], N: int) -> BinLinearMap:
    # !TODO! docstring for algorithm.generating_BinLinearMap_r
    return Blm(complete_basis_reversed(basis,N),N,N)


def generating_BinLinearMap(basis: list[BinWord], N: int) -> BinLinearMap:
    """Builds a full rank BinLinearMap such that the image of the inputs with the lowest MSBs corresponds to a specific basis.

    Args:
        basis: a list containing the integer representation of the binary vectors in the basis
        N: the intended input dimension for the output. It is assumed that the elements in basis are of length at most N.
    
    Returns
        A BinLinearMap L from N bits to N bits such that L(1 << i) = basis[i] for all i < len(basis) and such that L is a bijection.
    
    """
    return Blm(complete_basis(basis,N),N,N)


def BinLinearMap_from_masks(masks: list[BinWord], N: int, M: int) -> BinLinearMap:
    """Builds a BinLinearMap such that the image of the inputs with the lowest MSBs corresponds to a specific basis, and such that the other inputs are mapped to 0.

    Args:
        masks: a list containing the integer representation of the target binary vectors
        N: the intended input dimension for the output.
        M: the length of the output. It is assumed that the elements in basis are of length at most M.

    Returns:
        A BinLinearMap L from N bits to M bits such that L(1 << i) = basis[i] for all i < len(basis) and L(1<<i)=0 for all i >= len(basis).
    
    """
    return Blm(masks + [0 for _ in range(len(masks), N)], N, M)


def BinLinearMap_from_range_and_image(inputs: list[BinWord], outputs: list[BinWord], N:int, M: int) -> BinLinearMap:
    # !TODO! doctstring for algorithm.BinLinearMap_from_range_and_image
    return BinLinearMap_from_masks(outputs,M,N)*(generating_BinLinearMap(inputs,M).inverse())
    

# !SECTION! Solving Linear Systems


cdef class F2LinearSystem:
    """Implements a linear system of equation over F_2.

    Is optimized for the specific case of large systems where equations are added one by one.

    It has two modes, the selection being made during initialization.

    If `self.echelonize` is set to True, then only a set of independent equations is kept in memory, which can lead to a significant memory saving: instead of adding a redundant equation, it simply does nothing. However, it also means that the addition of new equations has a time complexity which is at worst linear in the total number of (independent) equations currently in the system. Since the rank is evaluated in real time (i.e., every time an equation is added), it could be used during the higher level computation itself.

    If `self.echelonize` is set to False, then equations are simply appended to the system without any further consideration. Redundant equations will not be simplified, but each addition is added in a constant (and small) time.
    
    Regardless of the value of `self.echelonize`, the solutions can only be obtained once the system is fully known.

    Under the hood, it uses the `cpp_F2LinearSystem` C++ class. Checks its documentation for more information.
    
    """

    def __init__(self, n_variables: int, echelonize: bool = False):
        """Initializes an F2LinearSystem instance.

        Args:
            n_variables(int): the number of variables that will intervene in the system
            echelonize(bool): whether the system should be echelonized on the fly
        
        """
        self.echelonize = echelonize
        self.cpp_ls = make_unique[cpp_F2LinearSystem](<int>n_variables, <bool>echelonize)

    
    def add_equation(self, variable_indices: list[BinWord]) -> bool:
        """Adds an equation to the system corresponding to the sum of the variables with the given indices. If the input is a list [0, 1, 3], then the equation x_0+x_1+x_3=0 is added to the system.

        If the system was initialized with `echelonize` set to True, then the running time is at worst proportional wit the rank of the current system as the system will get re-echelonized. If it was set to False, then the run time is constant.

        Args:
            variable_indices(list): a list containing the indices of the variables involved in the linear equation to be added to the system.

        Returns:
            If `echelonize` is False, then always returns True. Otherwise, returns True if the equation has increased the rank of the system, and False if the new equation was already in the span of the previous equations.
        
        """
        return dereference(self.cpp_ls).add_equation(<std_vector[BinWord]>variable_indices)

    
    def remove_solution(self, solution_indices: list[BinWord]) -> None:
        """Ensures that the space spanned by the Kernel of the system does not contain a specific value. If said value was indeed in the kernel, then the kernel dimension is decreased. If it actually was not in it, then nothing it will have no impact.

        The input describes the indices where the unwanted solution is set to 1. For example, if the input is [0, 1, 3], and if the system has 5 variables, then the solution (1,1,0,1,0) will be removed from the Kernel.

        In practice, this does not change the system of equations. Instead, a post-processing step is performed once the kernel is known to remove the (potential) contributions of all the solutions that have been removed using this method.

        Args:
           solution_indices(list): A list of indices corresponding to the positions set to 1 in the unwanted solution.
        
        """
        dereference(self.cpp_ls).remove_solution(<std_vector[BinWord]>solution_indices)

        
    def __len__(self) -> int:
        """Returns the number of equations currently in the system. If `self.echelonize`, then it is equal to the rank, otherwise, it directly corresponds to the number of equations that were added.

        Returns:
            An integer corresponding to the number of equations currently in the system.
        """
        return dereference(self.cpp_ls).size()

    
    def rank(self) -> int:
        """The current rank of the system, if it is known. If it isn't, throws an Exception.

        If you want to know the rank of a system that is not echelonized, you need to compute it by hand using the size of the kernel.
        
        Returns:
            If `self. echelonize` is set to True, then returns the current rank of the system. Otherwise, throws an Exception.
        """
        if self.echelonize:
            return dereference(self.cpp_ls).rank()
        else:
            raise Exception("Trying to get the rank of a non_echelonized system of equations")

        
    def kernel_as_bytes(self) -> list[bytes]:
        """Solves the system, removes unwanted solutions, and returns a basis of its kernel as list of `bytes` (each `bytes` being essentially a list of 8-bit unsigned char).

        Returns:
            A list where each entry is of type `bytes`. The bit with index i then corresponds to the value of x_i, and it can be obtained by computing e.g. `(y[i >> 3] >> (i & 0x7)) & 1`, where y is one of `bytes` contained in the output list.
        
        """
        return dereference(self.cpp_ls).kernel_as_bytes()


    def kernel_as_bits(self) -> list[bytes]:
        """Solves the system, removes unwanted solutions, and returns a basis of its kernel as list of `bytes` (each `bytes` being a list of that contains either 0 or 1.

        This method is less space-efficient than `kernel_as_bytes` because each bit needs a full byte of space; however, it means that the conversion between bits and bytes is done in the C++ world instead of the python world.

        Returns:
            A list where each entry is of type `bytes`. The bit with index i then corresponds to the value of x_i, and is simply equal to y[i], where y is one of `bytes` contained in the output list.
        
        """
        return dereference(self.cpp_ls).kernel_as_bits()
    
    
    def __str__(self) -> str:
        return dereference(self.cpp_ls).to_string().decode("UTF-8")
