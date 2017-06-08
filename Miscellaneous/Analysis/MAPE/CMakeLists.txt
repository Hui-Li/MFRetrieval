cmake_minimum_required(VERSION 2.8)
project(Coordinate)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -march=native -lgfortran -pthread")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++0x -march=native -lgfortran -pthread -I /home/hli2/libs/armadillo/usr/include")

#set(ARMADILLO_LIBRARY "/home/hli2/libs/armadillo/usr/lib/libarmadillo.so")
#set(ARMADILLO_INCLUDE_DIR "/home/hli2/libs/armadillo/usr/include")

find_package(Armadillo REQUIRED)

set(SOURCE_FILES main.cpp VectorElement.h Matrix.h)
add_executable(main ${SOURCE_FILES})
target_link_libraries(main ${ARMADILLO_LIBRARIES})