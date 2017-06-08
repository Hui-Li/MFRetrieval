cmake_minimum_required(VERSION 2.8)
project(Norm)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(SOURCE_FILES main.cpp Matrix.h VectorElement.h)
add_executable(Norm ${SOURCE_FILES})