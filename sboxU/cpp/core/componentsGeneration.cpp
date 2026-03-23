#include "./componentsGeneration.hpp"


cpp_S_box cpp_rand_invertible_S_box(cpp_PRNG & alea, unsigned int n)
{
    Lut perm = alea.get_permutation(1 << n);
    return cpp_S_box(perm);
}

cpp_S_box cpp_rand_S_box(
    cpp_PRNG & alea,
    unsigned int input_length,
    unsigned int output_length)
{
    Lut perm(1<<input_length, 0);
    for(BinWord x=0; x<perm.size(); x++)
        perm[x] = alea(0, 1 << output_length);
    return cpp_S_box(perm);
}
