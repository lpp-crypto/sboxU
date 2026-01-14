#include "../sboxU-cpp-module.hpp"
#include <cstdint>

void print_uint8_t(const uint8_t x)
{
    std::cout << " " << std::dec << std::setw(3) << (unsigned int)x;
}

void test_prng(cpp_PRNG & alea, uint8_t e)
{
    BinWord min = e/2, max = e;
    print_uint8_t(min);
    std::cout << " <= x < ";
    print_uint8_t(max);
    std::cout << "  : ";
    for(unsigned int i=0; i<20; i++)
        print_uint8_t(alea(min, max));
    std::cout << std::endl;
}
    

int main()
{
    std::cout << "Deterministic" << std::endl;
    for(uint8_t e=10; e<200; e+=40)
    {
        cpp_PRNG alea(std::vector<uint8_t>{{0, 2, e}});
        test_prng(alea, e);
    }
    std::cout << "\nRandomized seed" << std::endl;
    for(uint8_t e=10; e<200; e+=40)
    {
        cpp_PRNG alea(cpp_get_seed());
        test_prng(alea, e);
    }
    return 0;
}
