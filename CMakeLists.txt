# Time-stamp: <2017-10-06 19:09:15 lperrin>
#
# LICENSE

cmake_minimum_required(VERSION 3.2)


SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O2 -Wall -pedantic -std=c++11")


find_package(PythonLibs 2.7 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

ADD_LIBRARY(sboxu_cpp SHARED
  sboxu_cpp.cpp
  sboxu_cpp_utils.cpp
  sboxu_cpp_diff_lin.cpp
  sboxu_cpp_equiv.cpp)
TARGET_LINK_LIBRARIES(sboxu_cpp ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} pthread)
SET_TARGET_PROPERTIES(sboxu_cpp PROPERTIES PREFIX "")
