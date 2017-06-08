cmake_minimum_required(VERSION 2.8)
project(TimeAnalysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -march=native")

set(alg
    alg/BasicNormPrune.h)

set(util
    util/Base.h
    util/Monitor.h
    util/Calculator.h
    util/FileUtil.h)

set(structs
    structs/Matrix.h
    structs/ExtendMatrix.h
    structs/VectorElement.h
    structs/ExtendMatrixRow.h
    structs/IntMatrixRow.h)

add_library(lib SHARED ${alg} ${util} ${structs})
set_target_properties(lib PROPERTIES LINKER_LANGUAGE CXX)

set(SOURCE_FILES main.cpp)
add_executable(main ${SOURCE_FILES})
target_link_libraries(main lib)