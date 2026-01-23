cmake_minimum_required(VERSION 3.22)

project(
  sboxU_CPP
  VERSION 2.0
  LANGUAGES CXX)

SET(CMAKE_CXX_FLAGS "-O3 -march=native -std=c++20 -pthread -Wno-narrowing -fopenmp")



# !SECTION! Declaring the libraries

# !SUBSECTION! the CORE library

add_library(
  libcore
  ./core/anf.cpp
  ./core/s_box.cpp
  ./core/binLinearMap.cpp
  ./core/spectrum.cpp
  ./core/f2functions.cpp
  ./core/prng.cpp
  ./core/componentsGeneration.cpp
)

# !SUBSECTION! the STATISTICS library

add_library(
  libstatistics
  ./statistics/linear.cpp
  ./statistics/differential.cpp
  ./statistics/linear_structures.cpp
  ./statistics/boomerang.cpp
)
target_link_libraries(libstatistics PRIVATE libcore)


# !SUBSECTION! the ALGORITHMS library

add_library(
  libalgorithms
  ./algorithms/binLinearBasis.cpp
  ./algorithms/linearSystem.cpp
  ./algorithms/spaceSearch.cpp
)
target_link_libraries(libalgorithms PRIVATE libcore)


# !SUBSECTION! the CCZ library

add_library(
  libccz
  ./ccz/explore.cpp
  ./ccz/graph.cpp
  ./ccz/linear_equivalence.cpp
  ./ccz/zeroes.cpp
  ./ccz/linear_representative.cpp
  ./ccz/partition_preserving_linear_mapping/pplm.cpp
)
target_link_libraries(libccz PRIVATE libcore libstatistics libalgorithms)


# !SUBSECTION! the APN library

add_library(
  libapn
  ./apn/invariants.cpp
  ./apn/ccz_class.cpp
)
target_link_libraries(libapn PRIVATE libcore libstatistics libalgorithms libccz)



# !SUBSECTION! the main sboxU-CPP library

add_library(sboxU_CPP core statistics algorithms ccz) # these are folder names, not libary names!
Target_link_libraries(sboxU_CPP PRIVATE libcore libstatistics  libalgorithms libccz)


# !SECTION! Building tests

add_executable(test_f2Sbox ./tests/F2Sbox.cpp)
target_link_libraries(test_f2Sbox PRIVATE sboxU_CPP) 

add_executable(test_prng ./tests/prng.cpp)
target_link_libraries(test_prng PRIVATE sboxU_CPP) 
