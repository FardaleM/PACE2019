cmake_minimum_required (VERSION 2.6)
project (pace)

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=gnu++11 -g -march=native" )
add_library(vc main.cpp)
add_executable(submission main.cpp)
target_link_libraries (submission)

