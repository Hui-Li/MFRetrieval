cmake_minimum_required(VERSION 2.8)
project(SparsityMapping)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -std=c++11 -pthread -march=native -ftree-vectorize")

## look for boost
find_package(Boost COMPONENTS thread program_options system filesystem)
if(NOT ${Boost_FOUND})
    SET(BOOST_ROOT ~/local) # default
    SET(Boost_NO_SYSTEM_PATHS ON) # force to use own build
    find_package(Boost COMPONENTS thread program_options system filesystem)
endif(NOT ${Boost_FOUND})
if(Boost_FOUND)
    message(STATUS "Boost found")
    message(STATUS ${Boost_LIBRARIES})
    include_directories(${Boost_INCLUDE_DIR})
    link_libraries(${Boost_LIBRARIES})
endif(Boost_FOUND)

set(alg
    alg/TernaryBaseSchema.h
    alg/DaryBaseSchema.h)

set(util
    util/Base.h
    util/Monitor.h
    util/Calculator.h
    util/CompareUtil.h
    util/FileUtil.h)

set(structs
    structs/Matrix.h
    structs/VectorElement.h
    structs/FastHeap.h
    structs/Triplet.h)

add_library(SparsityMapping SHARED ${alg} ${util} ${structs})
set_target_properties(SparsityMapping PROPERTIES LINKER_LANGUAGE CXX)
link_libraries(SparsityMapping ${Boost_LIBRARIES})

add_executable(runNaive runNaive.cpp)
add_executable(runSparsityMapping runSparsityMapping.cpp)
add_executable(runSparsityMappingDary runSparsityMappingDary.cpp)