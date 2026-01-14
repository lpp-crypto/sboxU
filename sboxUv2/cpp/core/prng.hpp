#ifndef _PRNG_
#define _PRNG_

#include "../common.hpp"
#include "./f2functions.hpp"
#include <cstdint>
#include <random>


#define SEED_LENGTH 16

class cpp_PRNG
{
private:
    std::mt19937_64 gen;
    
public:
    /**
      Initializes a PRNG instance with the given string.

      The PRNG is completely deterministic: if the same string is used twice, then the sequence generated will be the same each time. By forcing an explicit seeding, we hope to foster reproducibility.

      
      @param seed A string that is used to initialize the instance.
     */
    inline cpp_PRNG(const std::vector<uint8_t> & seed)
    {
        if (seed.size() == 0)
            throw std::runtime_error("PRNG must be explicitely seeded");
        else 
        {
            std::seed_seq s(seed.begin(), seed.end());
            gen.seed(s);
        }
    }

    /** Returns a BinWord picked uniformly at random in [0, 2^64-1].
     */
    inline BinWord operator()()
    {
        return gen();
    }

           
    /** Returns a BinWord picked uniformly at random in [begin, end).

        The idea of the boundaries is to mimick the behaviour the python function `range`.

        @param min The minimum value that can be returned
        @param max The maximum value returned is max-1
     */
    inline BinWord operator()(BinWord begin, BinWord end)
    {
        if (end == begin)
            return begin;
        else if (end < begin)
            throw std::runtime_error("in PRNG(begin, end), end must be greater than begin");
        else
        {
            BinWord range = end - begin;
            if (cpp_hamming_weight(range) > 1)
            {
                BinWord
                    mask = (1 << (cpp_msb(range) + 2)) - 1,
                    result = range+1;
                while (result >= range)
                    result = gen() & mask;
                return begin + result;
            }
            else   // case where the range is a power of 2
            {
                return begin + (gen() & (range-1));
            }
        }
    }


    /** Returns a vector containing the integers in [0, ..., card) in a randomized order.

        Relies on a Fisher-Yates shuffles (see https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle).

        @param card The number of elements in the shuffled array.

        @return a vector of BinWord containing the shuffled integers.
     */
    inline Lut get_permutation(BinWord card)
    {
        Lut result(card, 0);
        for(unsigned int i=0; i<card; i++)
            result[i] = i;
        for(unsigned int i=0; i<card-1; i++)
        {
            unsigned int j = operator()(i, card);
            BinWord tmp = result[i];
            result[i] = result[j];
            result[j] = tmp;
        }
        return result;
    }
};

/** A helper function to obtain unique seeds, for example for tests.

    @return A vector of bytes of length 16 derived from the OS's default PRNG (on linux, /dev/urandom).
 */
inline std::vector<uint8_t> cpp_get_seed()
{
    std::vector<uint8_t> result(SEED_LENGTH, 0);
    std::random_device rd;
    std::uniform_int_distribution<int> dist(0, 255);
    std::cout << "seed used: {{";
    for(unsigned int i=0; i<SEED_LENGTH; i++)
    {
        result[i] = dist(rd);
        std::cout << std::dec << (unsigned int)result[i] << ",";
    }
    std::cout << "}" << std::endl;
    return result;
};
    
#endif
