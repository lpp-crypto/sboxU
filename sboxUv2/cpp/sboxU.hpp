#ifndef _SBOX_U_
#define _SBOX_U_


// !SECTION! Standard libraries

// STL containers
#include <vector>
#include <map>
#include <string>

// mathematics
#include <cmath>
#include <numeric>

// Input/Outputs
#include <iostream>
#include <sstream>
#include <iomanip>

// Exceptions
#include <exception>

// Multi-threading
#include <thread>

// Standard algorithms
#include <algorithm>


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

typedef std::vector<Integer> FpWord;

#endif
