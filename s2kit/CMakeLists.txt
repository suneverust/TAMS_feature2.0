# Afflication: TAMS, University of Hamburg
# E-Mail: bosun@informatik.uni-hamburg.de

# We gather the ESSENTIAL files for compute the spherical harmonics 
# in s2kit10, and write this CMakeLists.txt in order to let CMake to 
# configure and compile the project. 
 
# We use "FindFFTW3.cmake" to let CMake find the FFTW3_INCLUDE_DIR 
# and FFTW3_LIBRARIES automatically, please make sure 
# the "FindFFTW3.cmake" at the right place.


cmake_minimum_required (VERSION 2.8)

project(s2kit)

find_package (FFTW3 REQUIRED)

include_directories (${FFTW3_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/s2kit/include)

link_directories (${FFTW3_LIBRARIES})

add_library (s2kit 	makeweights.c
			primitive.c 
			naive_synthesis.c 
			pmls.c 
			cospmls.c 
			seminaive.c 
			FST_semi_memo.c)

target_link_libraries (s2kit ${FFTW3_LIBRARIES} -lm)
