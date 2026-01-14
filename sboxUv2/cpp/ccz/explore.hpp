#ifndef _CCZ_EXPLORE_
#define _CCZ_EXPLORE_

#include "../common.hpp"
#include "../core/include.hpp"
#include "./zeroes.hpp"
#include "./graph.hpp"


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

    @return A vector of S-boxes, each S-box being CCZ-equivalent s. Each EA-class in the CCZ-class of s will have at least one representative in this output.
 */
std::vector<cpp_S_box> cpp_enumerate_ea_classes(
    const cpp_S_box & s,
    const cpp_WalshZeroesSpaces &ws
    );

#endif
