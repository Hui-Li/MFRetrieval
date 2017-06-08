cmake_minimum_required(VERSION 2.8)
project(ValueDistribution)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -llapack -lblas")

find_package(Armadillo REQUIRED)

if(Armadillo_FOUND)
    message(STATUS "Armadillo found " ${ARMADILLO_LIBRARY_DIRS})
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    message(STATUS "Armadillo_LIBRARY_DIRS in  " ${ARMADILLO_LIBRARY_DIRS})
    message(STATUS "Armadillo_INCLUDE_DIRS in  " ${ARMADILLO_INCLUDE_DIRS})
endif()

set(SOURCE_FILES main.cpp Matrix.h VectorElement.h SVDDisAnalysis.h SVDDisAnalysis2.h DisAnalysis.h)
add_executable(main ${SOURCE_FILES})
target_link_libraries(main ${ARMADILLO_LIBRARIES})