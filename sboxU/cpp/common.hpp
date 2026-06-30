#ifndef _SBOX_U_
#define _SBOX_U_


// !SECTION! Standard libraries

#include <cstdint>

// STL containers
#include <vector>
#include <map>
#include <string>
#include <utility>

// mathematics
#include <cmath>
#include <numeric>

// Input/Outputs
#include <iostream>
#include <sstream>
#include <iomanip>

// Exceptions
#include <exception>
#include <stdexcept>


// Multi-threading
#include <thread>

// Standard algorithms
#include <algorithm>
#include <random>

// !SECTION! Special types 

// bytearray
typedef std::vector<uint8_t> Bytearray;

// signed
typedef int64_t Integer;
typedef std::vector<Integer> Row;

// unsigned
typedef uint64_t BinWord;
typedef std::vector<BinWord> Lut;

typedef std::vector< std::vector<std::vector<BinWord>>> Xtable;

typedef std::vector<uint64_t> FpWord;

#endif
