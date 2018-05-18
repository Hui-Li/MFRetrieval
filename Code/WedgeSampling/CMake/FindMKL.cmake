
# - Try to find mkl - the Intel Math Kernel Library
# Once done this will define
#  MKL_FOUND - System has mkl
#  MKL_INCLUDE_DIRS - The mkl include directories
#  MKL_LIBRARIES - The libraries needed to use mkl
#  MKL_DEFINITIONS - Compiler switches required for using mkl

if(NOT MKLROOT)
    set(MKLROOT /opt/intel/mkl)
endif()

find_path(MKL_INCLUDE_DIR
        NAMES mkl.h
        HINTS $ENV{MKLROOT}/include ${MKLROOT}/include
        PATHS)

find_library(MKL_LIBRARY
        NAMES mkl_rt
        HINTS $ENV{MKLROOT}/lib/intel64 ${MKLROOT}/lib/intel64 $ENV{MKLROOT}/lib ${MKLROOT}/lib
        PATHS)

set(MKL_LIBRARIES ${MKL_LIBRARY} )
set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(mkl  DEFAULT_MSG
        MKL_LIBRARY MKL_INCLUDE_DIR)

mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY)