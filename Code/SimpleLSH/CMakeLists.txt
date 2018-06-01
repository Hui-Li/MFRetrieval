cmake_minimum_required(VERSION 2.8)
project(SimpleLSH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -pthread -march=native -ftree-vectorize")

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


set(util
    util/Base.h
    util/Monitor.h
    util/Calculator.h
    util/Scanner.h
    util/Metric.h
    util/Topk.h
    util/MaxHeap.h)

set(structs
    structs/Matrix.h)

set(lsh
    lsh/psdLSH.h)

add_library(SimpleLSH SHARED ${util} ${structs} ${lsh})
set_target_properties(SimpleLSH PROPERTIES LINKER_LANGUAGE CXX)
link_libraries(SimpleLSH ${Boost_LIBRARIES})

add_executable(runSimpleLSH runSimpleLSH.cpp)