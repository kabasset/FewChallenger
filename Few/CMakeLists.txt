CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Few)

elements_depends_on_subdirs(ElementsKernel)

find_package(GSL)
find_package(HDF5)

elements_add_library(Few src/lib/*.cpp
                    INCLUDE_DIRS ElementsKernel GSL HDF5
                    LINK_LIBRARIES ElementsKernel GSL HDF5
                    PUBLIC_HEADERS Few)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                      EXECUTABLE Few_Amplitude_test
                      LINK_LIBRARIES Few
                      TYPE Boost)
