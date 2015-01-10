# Author: Bo Sun
# Afflication: TAMS, University of Hamburg
# E-Mail: bosun@informatik.uni-hamburg.de
#         user_mail@QQ.com

# We gather the ESSENTIAL files for compute the spherical harmonics 
# in s2kit10, and write this CMakeLists.txt in order to let CMake to 
# configure and compile the project. 
 
# We use "FindFFTW3.cmake" to let CMake find the FFTW3_INCLUDE_DIR 
# and FFTW3_LIBRARIES automatically, please make sure 
# the "FindFFTW3.cmake" at the right place.

# We use "FindGSL.cmake" to let CMake find the
# ${GSL_INCLUDE_DIRES} ${GSLCBLAS_INCLUDE_DIRS}
# and ${GSL_LIBRARIES} ${GSLBLAS_LIBRARIES},
# please make sure the "FindGSL.cmake" at the right place.

# Example: In my computer, the "FindFFTW3.cmake"
# and "FindGSL.cmake" are stored in '/usr/share/cmake-2.8/Modules'


cmake_minimum_required (VERSION 2.8)

project (tams_feature2.0)

find_package (FFTW3 REQUIRED)
include_directories (${FFTW3_INCLUDE_DIR})
link_directories (${FFTW3_LIBRARIES})

find_package(PCL 1.7 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

include_directories (${PROJECT_SOURCE_DIR}/s2kit/include)
add_subdirectory (${PROJECT_SOURCE_DIR}/s2kit)

add_library (tams_s2_semi_memo_for tams_s2_semi_memo_for.cpp)
target_link_libraries (tams_s2_semi_memo_for ${FFTW3_LIBRARIES} s2kit -lm)

add_executable (tams_feature2.0 tams_feature.cpp)

target_link_libraries (tams_feature2.0  ${FFTW3_LIBRARIES}
                                        ${PCL_LIBRARIES}
                                        s2kit
                                        tams_s2_semi_memo_for -lm)
