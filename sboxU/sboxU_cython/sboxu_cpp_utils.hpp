/* Time-stamp: <2021-08-11 15:08:37 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_UTILS_H_
#define _SBOXU_CPP_UTILS_H_

#include "sboxu_cpp.hpp"


// !SECTION! Basic functions 

// !SUBSECTION! Boolean operations 

/* @return the exclusive or of x and y. */ 
BinWord oplus_cpp(BinWord x, BinWord y) ;

/* @return the Hamming weight of x. */ 
unsigned int hamming_weight_cpp(BinWord x) ;

/* @return the parity of x, i.e. its hamming weight modulo 2. */ 
BinWord parity_cpp(BinWord x) ;

/* @return the scalar product of x and y. */ 
BinWord scal_prod_cpp(BinWord x, BinWord y) ;


/* @return the component $x \mapsto b . f(x)$. */ 
std::vector<BinWord> component_cpp(BinWord a, std::vector<BinWord> f);


// !SUBSECTION! Generating and testing permutations 

/* @return the LUT of a permutation of [0, 2^n) picked uniformly at
 * random using a PRNG from the C++ STL. */
Sbox random_permutation_cpp(unsigned int n);

/* @return true if and only if the LUT s contains all elements in [0,
 * s.size()). */
bool is_permutation_cpp(Sbox s);

/* @return the LUT of the functional inverse of s. */ 
Sbox inverse_cpp(Sbox s);


// !SECTION! Rank of sets of vectors

/* @return interpreting the elements in l as binary vectors of length
 * at most 64, returns the rank of this set of vectors.
 */
Integer rank_of_vector_set_cpp(std::vector<BinWord> l);


// !SECTION! Basic verifications

class NotSboxSized: public std::exception
{
public:
    unsigned int length;
    NotSboxSized(unsigned int l)
    {
        length = l;
    }
    
    const char* what() const throw()
    {
        std::string msg = "List has length ";
        msg += std::to_string(length);
        msg += " but length should be a power of 2";
        return msg.c_str();
    }
};

/* If the length of s is a power of 2, does nothing. Otherwise, throws
 * a NotSboxSized exception.
 */ 
void check_length_cpp(Sbox s);


#endif /* _SBOXU_CPP_UTILS_H_ */
