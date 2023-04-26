CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Challenger)

elements_depends_on_subdirs(ElementsKernel LinxCore LinxRun Few)

find_package(GSL)

elements_add_library(Challenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel LinxCore LinxRun Few GSL
                     LINK_LIBRARIES ElementsKernel LinxCore LinxRun Few GSL
                     PUBLIC_HEADERS Challenger)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                     EXECUTABLE Challenger_Amplitude_test
                     LINK_LIBRARIES Challenger
                     TYPE Boost)
elements_add_unit_test(Interpolant tests/src/Interpolant_test.cpp 
                     EXECUTABLE Challenger_Interpolant_test
                     LINK_LIBRARIES Challenger
                     TYPE Boost)
elements_add_unit_test(Spline tests/src/Spline_test.cpp 
                     EXECUTABLE Challenger_Spline_test
                     LINK_LIBRARIES Challenger
                     TYPE Boost)

elements_add_executable(Challenge src/program/Challenge.cpp
                     INCLUDE_DIRS ElementsKernel Challenger
                     LINK_LIBRARIES ElementsKernel Challenger)

elements_install_conf_files()
