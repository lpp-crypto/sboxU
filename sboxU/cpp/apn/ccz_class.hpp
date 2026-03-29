#ifndef _APN_CCZ_CLASS_
#define _APN_CCZ_CLASS_


#include "../common.hpp"
#include "../core/include.hpp"
#include "../ccz/include.hpp"
#include "./invariants.hpp"
#include "./ortho_derivative.hpp"


// Returns all EA automorphisms of the quadratic APN function s, using the
// ortho-derivative method.  Each entry is a pair (L, delta) where:
//
//   - L is the upper-triangular 2n×2n matrix encoding the automorphism,
//     with block decomposition  [L_B_T, L_C*L_A ; 0, L_A]:
//       block_A = L_B_T  — output transformation B of the EA automorphism
//       block_B = L_A    — inverse of input transformation A (so A = L_A^{-1})
//
//   - delta is the input shift: the automorphism satisfies
//       L_B_T(s(x)) = s(L_A^{-1}(x) XOR delta) XOR L_C(x) XOR C_0
//     for constants L_C (linear) and C_0.
//
// To compare with the lower-triangular (r, a) pairs returned by
// ea_equivalences / cpp_equivalences_from_lat, use the canonical
// (A_affine, B) form:
//   A_affine(x) = block_B(L)^{-1}(x) XOR delta  (vs. block_A(r)(x) XOR a)
//   B(x)        = block_A(L)(x)                  (vs. block_B(r)(x))
//
// Note: the lower-triangular convention has block_C = 0, block_D != 0;
//       the upper-triangular convention here has block_C != 0, block_D = 0.
std::vector<std::pair<cpp_F2AffineMap, BinWord>> cpp_automorphisms_from_ortho_derivative(
    const cpp_S_box & s,
    const unsigned int n_threads
    );

std::vector<cpp_F2AffineMap> cpp_ea_mappings_from_ortho_derivative(
    const cpp_S_box & s,
    const cpp_S_box & s_prime,
    const unsigned int n_threads
    );
    
std::vector<cpp_S_box> cpp_enumerate_ea_classes_quadratic_apn(
    const cpp_S_box &s,
    const unsigned int n_threads
    );


cpp_S_box cpp_ccz_equivalent_quadratic_function(
    const cpp_S_box & s,
    const unsigned int n_threads
    );


#endif
