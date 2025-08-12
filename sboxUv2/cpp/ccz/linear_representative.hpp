#ifndef _CCZ_LINEAR_EQUIV_
#define _CCZ_LINEAR_EQUIV_

#include "../sboxU.hpp"
#include "../s_box.hpp"


typedef std::pair<uint64_t, uint64_t> IOpair;


// !SECTION! Python-facing functions

// std::vector<Sbox> cpp_linear_equivalence(const Sbox f, const Sbox g, bool all_mappings);

cpp_S_box cpp_class_representative(const cpp_S_box & f);

    
// !SECTION! Exceptions needed to implement guess and determine 


// This class is used to track contradiction in guess & determine
// algorithms.
class ContradictionFound: public std::exception
{
public:
    std::string msg;
    ContradictionFound(BinWord x, BinWord y)
    {
        msg = ""; // let's not waste time initializing a string
    }
    
    const char* what() const throw()
    {
        return msg.c_str();
    }
};



// !SECTION! Helper classes 
// !SUBSECTION! Classes for guesses of linear permutations 

// This class stores a guess for a linear permutation as used in the
// algorithm checking linear equivalence.
class LEguess
{
private:
    unsigned int target_size;
    std::map<BinWord, BinWord> partial_lut;
    std::vector<IOpair> latest_entries;
    unsigned int min_unset;
public:
    LEguess(const unsigned int _target_size);
    ~LEguess();
    std::map<BinWord, bool> is_set;
    std::vector<IOpair> add_entry(const IOpair e);
    bool is_entry_set(const BinWord x);
    unsigned int min_u();
    bool complete();
    BinWord img(const BinWord x);
    Lut lut();
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
    LEguessIterator(
        const unsigned int _target_size,
        const std::vector<IOpair> _constraints
        );
    bool prepare_successfully();
    LEguess get_prepared_guess();
    LEguessIterator deeper_guess_gen();
};


// !SUBSECTION! Class to store the guess of the representative of a non-linear permutation in its linear equivalence class.

class LErepr
{
private:
    cpp_S_box f;
    cpp_S_box f_inv;
    std::map<BinWord, BinWord> partial_lut;
    std::map<BinWord, BinWord> partial_inv_lut;
    std::map<BinWord, bool> is_set;
    std::map<BinWord, bool> is_inv_set;
    unsigned int contiguous_length;
    Lut best_guess;
    LEguess a;
    LEguess b;
    LEguess a_inv;
    LEguess b_inv;
public:
    LErepr(const cpp_S_box _f);
    void add_entry(IOpair e);
    bool current_is_greater_than_best();
    unsigned int min_n_a();
    unsigned int min_n_b();
    unsigned int min_u_a();
    unsigned int min_u_b();
    bool update_state();
    void initialize();
    Lut lut();
};



#endif /* _SBOXU_CPP_EQUIV_H_ */
