cmake_minimum_required(VERSION 2.8)
project(BoundAnalysis)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -std=c++0x")

set(alg
    alg/int/MapIntBasicPrune.h)

set(util
    util/Base.h
    util/Calculator.h)

set(structs
    structs/Matrix.h
    structs/ExtendMatrix.h
    structs/VectorElement.h
    structs/ExtendMatrixRow.h
    structs/IntMatrixRow.h)

add_library(lib SHARED ${alg} ${util} ${structs})
set_target_properties(lib PROPERTIES LINKER_LANGUAGE CXX)
set(SOURCE_FILES main.cpp IntScale.h)
add_executable(main ${SOURCE_FILES})
target_link_libraries(main lib)