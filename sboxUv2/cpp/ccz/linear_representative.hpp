#ifndef _CCZ_LINEAR_EQUIV_
#define _CCZ_LINEAR_EQUIV_

#include "../sboxU.hpp"
#include "../core/include.hpp"




std::tuple<cpp_S_box, cpp_BinLinearMap, cpp_BinLinearMap> cpp_class_representative(const cpp_S_box & f);


// the following class was used by the old linear representative implementation
// typedef std::pair<uint64_t, uint64_t> IOpair;

// class LErepr
// {
// private:
//     cpp_S_box f;
//     cpp_S_box f_inv;
//     std::map<BinWord, BinWord> partial_lut;
//     std::map<BinWord, BinWord> partial_inv_lut;
//     std::map<BinWord, bool> is_set;
//     std::map<BinWord, bool> is_inv_set;
//     unsigned int contiguous_length;
//     Lut best_guess;
//     LEguess a;
//     LEguess b;
//     LEguess a_inv;
//     LEguess b_inv;
// public:
//     LErepr(const cpp_S_box _f);
//     void add_entry(IOpair e);
//     bool current_is_greater_than_best();
//     unsigned int min_n_a();
//     unsigned int min_n_b();
//     unsigned int min_u_a();
//     unsigned int min_u_b();
//     bool update_state();
//     void initialize();
//     Lut lut();
// };




#endif /* _SBOXU_CPP_EQUIV_H_ */
