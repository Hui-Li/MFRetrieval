cmake_minimum_required(VERSION 3.5)
project(PositiveDist)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")

set(SOURCE_FILES main.cpp Matrix.h VectorElement.h)
add_executable(PositiveDist ${SOURCE_FILES})