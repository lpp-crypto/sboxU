#ifndef _SBOX_U_
#define _SBOX_U_


// !SECTION! Standard libraries

// STL containers
#include <vector>
#include <map>
#include <string>
#include <tuple>

// Input/Outputs
#include <iostream>
#include <iomanip>

// Exceptions
#include <exception>

// Multi-threading
#include <thread>


// !SECTION! Special types 


typedef int64_t Integer;
typedef std::vector<Integer> Row;
typedef uint64_t BinWord;




// !SECTION! Components

#include "f2functions.hpp"
#include "s_box.hpp"
#include "statistics/spectrum.hpp"

#endif
