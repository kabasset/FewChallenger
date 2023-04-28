CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(MemoryChallenger)

elements_depends_on_subdirs(ElementsKernel LinxCore LinxRun)

find_package(GSL)

elements_add_library(MemoryChallenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel LinxCore LinxRun GSL
                     LINK_LIBRARIES ElementsKernel LinxCore LinxRun GSL
                     PUBLIC_HEADERS MemoryChallenger)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                     EXECUTABLE Challenger_Amplitude_test
                     LINK_LIBRARIES MemoryChallenger
                     TYPE Boost)
elements_add_unit_test(Interpolant tests/src/Interpolant_test.cpp 
                     EXECUTABLE Challenger_Interpolant_test
                     LINK_LIBRARIES MemoryChallenger
                     TYPE Boost)
