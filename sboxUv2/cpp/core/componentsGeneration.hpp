#ifndef _COMPONENTS_GENERATION_
#define _COMPONENTS_GENERATION_

#include "../common.hpp"
#include "./s_box.hpp"
#include "./f2functions.hpp"
#include "./binLinearMap.hpp"
#include "./prng.hpp"


/** Returns a cpp_S_box instance corresponding to an invertible transformation picked uniformly at random.

    @param alea The source of randomness to use.
    @param n The bit length of both the input and output of the transformation.

    @return a cpp_S_box instance whose LUT contains each integer in {0, ..., 2^n-1} exactly once.
 */
cpp_S_box cpp_rand_invertible_S_box(cpp_PRNG & alea, unsigned int n);


/** Returns a cpp_S_box instance corresponding to a function whose outputs are picked independently and uniformly at random.

    @param alea The source of randomness to use.
    @param input_length The bit length of the input.
    @param output_length The bit length of the output.

    @return a cpp_S_box instance whose LUT contains 2^input_length integers of {0, ..., 2^output_length-1}, each picked uniformly at random (independently from each other).
 */
cpp_S_box cpp_rand_S_box(
    cpp_PRNG & alea,
    unsigned int input_length,
    unsigned int output_length);


#endif
