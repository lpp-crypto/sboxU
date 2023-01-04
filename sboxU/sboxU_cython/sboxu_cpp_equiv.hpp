/* Time-stamp: <2023-01-04 15:52:05 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_EQUIV_H_
#define _SBOXU_CPP_EQUIV_H_

#include "sboxu_cpp.hpp"



// !SECTION! Python-facing functions

std::vector<Sbox> linear_equivalence_cpp(const Sbox f, const Sbox g, bool all_mappings);
Sbox le_class_representative_cpp(const Sbox f);

    
// !SECTION! Exceptions needed to implement guess and determine 


// This class is used to track contradiction in guess & determine
// algorithms.
class ContradictionFound: public std::exception
{
public:
    std::string msg;
    ContradictionFound(uint64_t x, uint64_t y)
    {
        msg = "Contradiction found: " + std::to_string(x) + " -> " + std::to_string(y) ;
    }
    
    const char* what() const throw()
    {
        return msg.c_str();
    }
};


// !SECTION! Classes for guesses of linear permutations 

// This class stores a guess for a linear permutation as used in the
// algorithm checking linear equivalence.
class LEguess
{
private:
    unsigned int target_size;
    std::map<uint64_t, uint64_t> partial_lut;
    std::vector<IOpair> latest_entries;
    unsigned int min_unset;
public:
    LEguess(const unsigned int _target_size);
    ~LEguess();
    std::map<uint64_t, bool> is_set;
    std::vector<IOpair> add_entry(const IOpair e);
    bool is_entry_set(const uint64_t x);
    unsigned int min_u();
    bool complete();
    uint64_t img(const uint64_t x);
    Sbox lut();
};

class LEguessIterator
{
private:
    unsigned int target_size;
    unsigned int guess_bit_length;
    unsigned long int x;
    unsigned long int y;
    unsigned long int guess_mask;
    std::vector<IOpair> constraints;
    LEguess base_guess;
    LEguess prepared_guess;
public:
    LEguessIterator(const unsigned int _target_size,
                    const std::vector<IOpair> _constraints);
    bool prepare_successfully();
    LEguess get_prepared_guess();
    LEguessIterator deeper_guess_gen();
};


// !SECTION! Class to store the guess of the representative of a
// !non-linear permutation in its linear equivalence class.

class LErepr
{
private:
    Sbox f;
    Sbox f_inv;
    std::map<uint64_t, uint64_t> partial_lut;
    std::map<uint64_t, uint64_t> partial_inv_lut;
    std::map<uint64_t, bool> is_set;
    std::map<uint64_t, bool> is_inv_set;
    unsigned int contiguous_length;
    Sbox best_guess;
    LEguess a;
    LEguess b;
    LEguess a_inv;
    LEguess b_inv;
public:
    LErepr(const Sbox _f);
    void add_entry(IOpair e);
    bool current_is_greater_than_best();
    unsigned int min_n_a();
    unsigned int min_n_b();
    unsigned int min_u_a();
    unsigned int min_u_b();
    bool update_state();
    void initialize();
    Sbox lut();
};



#endif /* _SBOXU_CPP_EQUIV_H_ */
