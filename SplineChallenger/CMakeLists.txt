CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(SplineChallenger)

elements_depends_on_subdirs(ElementsKernel LinxCore Splider MemoryChallenger)

elements_add_library(SplineChallenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel LinxCore Splider MemoryChallenger
                     LINK_LIBRARIES ElementsKernel LinxCore Splider MemoryChallenger
                     PUBLIC_HEADERS SplineChallenger)

elements_add_unit_test(Amplitude tests/src/Amplitude_test.cpp 
                    EXECUTABLE SplineChallenger_Amplitude_test
                    LINK_LIBRARIES SplineChallenger
                    TYPE Boost)
