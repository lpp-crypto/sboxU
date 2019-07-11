/* Time-stamp: <2019-07-10 17:12:49 lperrin>
 *
 * LICENSE
 */ 

#ifndef _SBOXU_CPP_H_
#define _SBOXU_CPP_H_


// !SECTION! Includes

// STL containers
#include <vector>
#include <utility>
#include <map>
#include <string>

// Basic libraries
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <random>
#include <exception>
#include <stdexcept>
#include <iostream>

// Using multi-threading
#include <thread>


// Seting up for boost.python to wrap the following functions
 #include <boost/python.hpp>
using namespace boost::python;



// !SECTION! Special types 

typedef int32_t Integer;
typedef uint32_t BinWord;
typedef std::vector<BinWord> Sbox;
typedef std::pair<uint32_t, uint32_t> IOpair;


// !SECTION! Including sub-modules

#include "sboxu_cpp_utils.hpp"
#include "sboxu_cpp_diff_lin.hpp"
#include "sboxu_cpp_equiv.hpp"
#include "sboxu_cpp_ccz.hpp"




#endif /* _SBOXU_CPP_H_ */
