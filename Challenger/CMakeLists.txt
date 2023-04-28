CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(Challenger)

elements_depends_on_subdirs(ElementsKernel LinxRun Few MemoryChallenger SplineChallenger)

find_package(GSL)

elements_add_library(Challenger src/lib/*.cpp
                     INCLUDE_DIRS ElementsKernel LinxRun Few MemoryChallenger SplineChallenger
                     LINK_LIBRARIES ElementsKernel LinxRun Few MemoryChallenger SplineChallenger
                     PUBLIC_HEADERS Challenger)

elements_add_executable(Challenge src/program/Challenge.cpp
                     INCLUDE_DIRS ElementsKernel Challenger
                     LINK_LIBRARIES ElementsKernel Challenger)

elements_install_conf_files()
