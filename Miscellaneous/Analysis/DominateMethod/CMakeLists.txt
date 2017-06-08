cmake_minimum_required(VERSION 3.5)
project(DominateMethod)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x")

set(SOURCE_FILES main.cpp Matrix.h VectorElement.h Calculator.h)
add_executable(DominateMethod ${SOURCE_FILES})