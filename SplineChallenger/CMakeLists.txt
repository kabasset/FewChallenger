CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(SplineChallenger)

elements_depends_on_subdirs(ElementsKernel LinxCore)

find_package(GSL)

elements_add_library(SplineChallenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel LinxCore GSL
                     LINK_LIBRARIES ElementsKernel LinxCore GSL
                     PUBLIC_HEADERS SplineChallenger)

elements_add_unit_test(Interpolant tests/src/Interpolant_test.cpp 
                    EXECUTABLE SplineChallenger_Interpolant_test
                    LINK_LIBRARIES SplineChallenger
                    TYPE Boost)

elements_add_unit_test(Spline tests/src/Spline_test.cpp 
                    EXECUTABLE SplineChallenger_Spline_test
                    LINK_LIBRARIES SplineChallenger
                    TYPE Boost)
