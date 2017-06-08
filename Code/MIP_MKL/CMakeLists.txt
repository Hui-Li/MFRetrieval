cmake_minimum_required(VERSION 2.8)
project(MIP_BatchNaive)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -std=c++11 -march=native -ftree-vectorize -lgfortran -pthread")

set(alg
    alg/ArmaNaive.h
    alg/MKLNaive.h)

set(structs
    structs/VectorElement.h)

set(util
    util/Base.h
    util/Monitor.h
    util/FileUtil.h
    util/Conf.h
    util/Logger.h)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMake/")

find_package(MKL REQUIRED)

if(MKL_FOUND)
    message(STATUS "MKL found " ${MKL_INCLUDE_DIRS})
    include_directories(${MKL_INCLUDE_DIRS})
    link_directories(${MKL_LIBRARIES})
endif()

find_package(Armadillo REQUIRED)

if(Armadillo_FOUND)
    message(STATUS "Armadillo found " ${ARMADILLO_LIBRARY_DIRS})
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    message(STATUS "Armadillo_LIBRARY_DIRS in  " ${ARMADILLO_LIBRARY_DIRS})
    message(STATUS "Armadillo_INCLUDE_DIRS in  " ${ARMADILLO_INCLUDE_DIRS})
endif()

find_package(Boost COMPONENTS serialization program_options)

if(NOT ${Boost_FOUND})
    SET(BOOST_ROOT ~/local) # default
    SET(Boost_NO_SYSTEM_PATHS ON) # force to use own build
    find_package(Boost COMPONENTS serialization program_options REQUIRED)
endif(NOT ${Boost_FOUND})

if(Boost_FOUND)
    message(STATUS "Boost found " ${Boost_LIBRARY_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    include_directories(${Boost_INCLUDE_DIRS})
    message(STATUS "Boost_LIBRARY_DIRS in  " ${Boost_LIBRARY_DIRS})
    message(STATUS "Boost_INCLUDE_DIRS in  " ${Boost_INCLUDE_DIRS})
endif()

add_library(lib SHARED ${alg} ${util} ${structs})
set_target_properties(lib PROPERTIES LINKER_LANGUAGE CXX)

add_executable(main main.cpp)
target_link_libraries(main lib)
target_link_libraries(main ${MKL_LIBRARIES})
target_link_libraries(main ${ARMADILLO_LIBRARIES})
target_link_libraries(main ${Boost_LIBRARIES})