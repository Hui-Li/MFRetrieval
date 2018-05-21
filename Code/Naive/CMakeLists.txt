cmake_minimum_required(VERSION 2.8)
project(Naive)

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


set(structs
    structs/FastHeap.h
    structs/VectorElement.h)

set(util
    util/Monitor.h
    util/Base.h
    util/FileUtil.h)

add_library(Naive SHARED ${util} ${structs})
set_target_properties(Naive PROPERTIES LINKER_LANGUAGE CXX)
link_libraries(Naive ${Boost_LIBRARIES} ${MKL_LIBRARIES})

add_executable(runNaive runNaive.cpp)