#ifndef _CCZ_LINEAR_EQUIVALENCE_
#define _CCZ_LINEAR_EQUIVALENCE_

#include "./linear_representative.hpp"


std::vector<cpp_F2AffineMap> cpp_linear_equivalence_permutations(
    const cpp_S_box f,
    const cpp_S_box g,
    bool all_mappings);
    

typedef std::pair<BinWord, BinWord> IOpair;


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


#endif
