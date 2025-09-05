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
#include <iomanip>

// Exceptions
#include <exception>

// Multi-threading
#include <thread>

// Standard algorithms
#include <algorithm>


// !SECTION! Special types 

// signed
typedef int64_t Integer;
typedef std::vector<Integer> Row;
// unsigned
typedef uint64_t BinWord;
typedef std::vector<BinWord> Lut;

typedef std::vector<Integer> FpWord;

#endif
