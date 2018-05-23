/* Time-stamp: <2018-05-16 14:40:12 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_UTILS_H_
#define _SBOXU_CPP_UTILS_H_

#include "sboxu_cpp.hpp"
using namespace boost::python;


// !SECTION! Basic functions 

// !SUBSECTION! Boolean operations 

/* @return the exclusive or of x and y. */ 
BinWord oplus_cpp(BinWord x, BinWord y) ;

/* @return the Hamming weight of x. */ 
unsigned int hamming_weight(BinWord x) ;

/* @return the parity of x, i.e. its hamming weight modulo 2. */ 
BinWord parity(BinWord x) ;

/* @return the scalar product of x and y. */ 
BinWord scal_prod(BinWord x, BinWord y) ;


/* @return the component $x \mapsto b . f(x)$. */ 
std::vector<BinWord> component(BinWord a, std::vector<BinWord> f);


// !SUBSECTION! Generating and testing permutations 

/* @return the LUT of a permutation of [0, 2^n) picked uniformly at
 * random using a PRNG from the C++ STL. */
Sbox random_permutation_cpp(unsigned int n);

/* @return true if and only if the LUT s contains all elements in [0,
 * s.size()). */
bool is_permutation_cpp(Sbox s);

/* @return the LUT of the functional inverse of s. */ 
Sbox inverse(Sbox s);


// !SUBSECTION! C++/Python converstion 

/* @return a C++ vector of Integers containing the signed numbers in l.
 *
 * @exception throws an exception if an element in the list is not an integer.
 */
std::vector<Integer> lst_2_vec_Integer(const list& l);

/* @return a C++ vector of BinWord containing the binary representations of the unsigned integers in l.
 *
 * @exception throws an exception if an element in the list is not an integer.
 */
std::vector<BinWord> lst_2_vec_BinWord(const list& l);


/* @return a python list of Integers containing the signed numbers in l.
 */
list vec_2_lst_Integer(const std::vector<Integer> l);

/* @return a python list of Integers containing the numbers whose binary expansion are the words in l.
 */
list vec_2_lst_BinWord(const std::vector<BinWord> l);



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
void check_length(Sbox s);




#endif /* _SBOXU_CPP_UTILS_H_ */
