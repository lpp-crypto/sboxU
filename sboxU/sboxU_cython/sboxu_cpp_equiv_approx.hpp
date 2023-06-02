/* Time-stamp: <2023-01-04 15:52:05 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_EQUIV_APPROX_H_
#define _SBOXU_CPP_EQUIV_APPROX_H_

#include "sboxu_cpp.hpp"



// !SECTION! Python-facing functions

std::vector<Sbox> linear_equivalence_approx_cpp(const Sbox f,
						const Sbox g,
						const bool all_mappings,
						const unsigned int max_contradictions);

    
// !SECTION! Exceptions needed to implement guess and determine 


// This class is used to track contradiction in guess & determine
// algorithms, and handles a maximum number of contradiction
class Opponent
{
private:
    unsigned int max_contradictions;
    unsigned int n_contradictions;

public:
    Opponent(unsigned int _max_contradictions) :
        max_contradictions(_max_contradictions),
	n_contradictions(0)
    {};

    void reset()
    { 
        n_contradictions = 0;
    };
  
    void process_contradiction()
    {
        n_contradictions ++;
	if (n_contradictions > max_contradictions)
            throw ContradictionFound(0, 0);
    };
};

  
// !SECTION! Classes for guesses of linear permutations 

// This class stores a guess for a linear permutation as used in the
// algorithm checking linear equivalence.
class LEguessApprox
{
private:
    unsigned int target_size;
    Opponent opp;
    std::map<uint64_t, uint64_t> partial_lut;
    std::vector<IOpair> latest_entries;
    unsigned int min_unset;
public:
    LEguessApprox(const unsigned int _target_size, Opponent _opp);
    ~LEguessApprox();
    std::map<uint64_t, bool> is_set;
    std::vector<IOpair> add_entry(const IOpair e);
    bool is_entry_set(const uint64_t x);
    unsigned int min_u();
    bool complete();
    uint64_t img(const uint64_t x);
    Sbox lut();
};


class LEguessIteratorApprox
{
private:
    unsigned int target_size;
    unsigned int guess_bit_length;
    unsigned long int x;
    unsigned long int y;
    unsigned long int guess_mask;
    std::vector<IOpair> constraints;
    LEguessApprox base_guess;
    LEguessApprox prepared_guess;
    Opponent opp;
public:
    LEguessIteratorApprox(const unsigned int _target_size,
			  const std::vector<IOpair> _constraints,
			  Opponent _opp);
    bool prepare_successfully();
    LEguessApprox get_prepared_guess();
    LEguessIteratorApprox deeper_guess_gen();
};





#endif /* _SBOXU_CPP_EQUIV_APPROX_H_ */
