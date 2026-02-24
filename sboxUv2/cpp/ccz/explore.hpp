#ifndef _CCZ_EXPLORE_
#define _CCZ_EXPLORE_

#include "../common.hpp"
#include "../core/include.hpp"
#include "./zeroes.hpp"
#include "./graph.hpp"



// !SECTION! Basic functions
//
// ! The following functions simply operate on the graph of a given
// ! function in order either to obtain a new function by applying a
// ! linear permutation to their graph, or by first deriving these
// ! functions from the Walsh zeroes.
// 
// ! Understanding the content of this section of the code requires
// ! a good familiarity with [1], in particular reagarding the
// ! definition and importance of the Walsh zeroes.
//
// ! [1] Canteaut, A., & Perrin, L. (2019). On CCZ-equivalence,
// ! extended-affine equivalence, and function twisting. Finite Fields
// ! and Their Applications, 56, 209-246.


/** Applies the given CCZ mapping to the graph of the given function and returns the function thus obtained.

   If the mapping is not actually admissible for s, returns the empty S-box. 

    @param s The cpp_S_box representing the function analysed.
    @param L The cpp_BinLinearMap to apply to the graph of s

    @return The function with graph L(G_s), where G_s is the graph of s. If L is not admissible for s, then returns the empty S-box.
 */
cpp_S_box cpp_ccz_equivalent_function(
    const cpp_S_box & s,
    const cpp_BinLinearMap & L
    );


/** Iterates over representatives of the extended-affine equivalence classes in the CCZ class of the given S-box.

    Every EA-classes in the CCZ-class of s will have at least one representative in the output of this function. However, there is a priori no guarantee that this representative will be unique: each EA-class is likely to be represented multiple times.

    It works by first generating the Walsh zeroes of s, then finding the vector spaces of dimension n in it, building the corresponding CCZ linear permutation, and finally using the function cpp_ccz_equivalent_function to apply this mapping to the graph of s.

    @param s The cpp_S_box representing the function analysed.
    @param n_threads The number of threads to use when searching for the vector spaces.

    @return A vector of S-boxes, each S-box being CCZ-equivalent s. Each EA-class in the CCZ-class of s will have at least one representative in this output.
 */
std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const unsigned int n_threads
    );


/** Iterates over representatives of the extended-affine equivalence classes in the CCZ class of the given S-box.

    Every EA-classes in the CCZ-class of s will have at least one representative in the output of this function. However, there is a priori no guarantee that this representative will be unique: each EA-class is likely to be represented multiple times.

    Unlike cpp_enumerate_ea_classes(const cpp_S_box & s, const unsigned int n_threads), this function takes the spaces in the Walsh zeroes as an input, thus skipping its first steps. It thus starts by building the corresponding CCZ linear permutations, and then uses the function cpp_ccz_equivalent_function to apply each of these mappings to the graph of s.

    @param s The cpp_S_box representing the function analysed.
    @param ws The spaces contained in the Walsh zeroes of s. No verifications are made: if the object supplied doesn't actually corresponds to s, it may result in a crash---or in any other kind of undefined behaviours.

    @return A vector of S-boxes, each S-box being CCZ-equivalent to s. Each EA-class in the CCZ-class of s will have at least one representative in this output.
 */
std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const cpp_WalshZeroesSpaces &ws
    );


// !SECTION! CCZ-equivalence to a permutation

// ! Using again the results in [1], it is possible to specifically
// ! target the permutations in the CCZ-class of a function. The
// ! functions below are doing just that, for the most general form of
// ! CCZ-equivalence, and for the specific case of extended-affine
// ! equivalence as well.

// !SUBSECTION! The general case

/** Finds at least one representative of each affine equivalence class of
   permutations contained in the CCZ-class of the given S-box.

   It works by identifying pairs U,U' of spaces in the Walsh zeroes of the
   given S-box such that their direct sum is the full space, then creating a
   mapping L such that L(V) = U and L(V^T) = U' (using the notations of [1]),
   and finally applying L^T to the graph of the given S-box.

   [1] Canteaut, A., & Perrin, L. (2019). On CCZ-equivalence, extended-affine
   equivalence, and function twisting. Finite Fields and Their Applications, 56,
   209-246.

   @param s The cpp_S_box representing the function analysed.
   @param ws The spaces contained in the Walsh zeroes of s. No verifications
   are made: if the object supplied doesn't actually corresponds to s, it may
   result in a crash---or in any other kind of undefined behaviours.

   @return A vector of S-boxes, each S-box being CCZ-equivalent to s. Each
   EA-class in the CCZ-class of s will have at least one representative in this
   output.
*/
std::vector<cpp_S_box> cpp_enumerate_permutations_in_ccz_class(
    const cpp_S_box &s,
    const cpp_WalshZeroesSpaces &ws
    );


/** Finds at least one representative of each affine equivalence class of
   permutations contained in the CCZ-class of the given S-box.

   It works by first computing the Walsh zeroes of the function, and then calling the other cpp_enumerate_permutations_in_ccz_class function on the result.

   @param s The cpp_S_box representing the function analysed.
   @param n_threads The number of threads to use when searching for the vector spaces.

   @return A vector of S-boxes, each S-box being CCZ-equivalent to s. Each
   EA-class in the CCZ-class of s will have at least one representative in this
   output.
 */
std::vector<cpp_S_box> cpp_enumerate_permutations_in_ccz_class(
    const cpp_S_box & s,
    const unsigned int n_threads
    );


// !SUBSECTION! The particular case of EA-equivalence


// !TODO! write the skeleton of functions for EA-equivalence to a permutation 

#endif
