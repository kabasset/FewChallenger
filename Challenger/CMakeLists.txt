CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Challenger)

elements_depends_on_subdirs(ElementsKernel LinxCore)

find_package(GSL)

elements_add_library(Challenger src/lib/*.cpp
                    INCLUDE_DIRS ElementsKernel LinxCore GSL
                    LINK_LIBRARIES ElementsKernel LinxCore GSL
                    PUBLIC_HEADERS Challenger)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                    EXECUTABLE Challenger_Amplitude_test
                    LINK_LIBRARIES Challenger
                    TYPE Boost)
