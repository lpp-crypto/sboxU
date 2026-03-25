#include "./prng.hpp"

Lut
cpp_PRNG::get_permutation(BinWord card)
{
  Lut result(card, 0);
  for (unsigned int i = 0; i < card; i++) {
    result[i] = i;
  }
  for (unsigned int i = 0; i < card - 1; i++) {
    unsigned int j = operator()(i, card);
    BinWord tmp = result[i];
    result[i] = result[j];
    result[j] = tmp;
  }
  return result;
}

static void
cpp_print_seed(const Bytearray& seed)
{
  std::cout << "seed used: {{";
  for (unsigned int i = 0; i < seed.size(); i++) {
    std::cout << std::dec << (unsigned int)seed[i] << ",";
  }
  std::cout << "}" << std::endl;
}

Bytearray
cpp_get_seed()
{
  Bytearray result(SEED_LENGTH, 0);
  std::random_device rd;
  std::uniform_int_distribution<int> dist(0, 255);
  for (unsigned int i = 0; i < SEED_LENGTH; i++) {
    result[i] = dist(rd);
  }
  cpp_print_seed(result);
  return result;
}
