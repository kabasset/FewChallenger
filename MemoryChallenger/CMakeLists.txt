CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(MemoryChallenger)

elements_depends_on_subdirs(ElementsKernel Linx)

find_package(GSL)

elements_add_library(MemoryChallenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel Linx GSL
                     LINK_LIBRARIES ElementsKernel Linx GSL
                     PUBLIC_HEADERS MemoryChallenger)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                     EXECUTABLE MemoryChallenger_Amplitude_test
                     LINK_LIBRARIES MemoryChallenger
                     TYPE Boost)
elements_add_unit_test(Interpolant tests/src/Interpolant_test.cpp 
                     EXECUTABLE MemoryChallenger_Interpolant_test
                     LINK_LIBRARIES MemoryChallenger
                     TYPE Boost)
