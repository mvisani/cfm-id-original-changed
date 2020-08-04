##########################################################################
# Mass Spec Prediction and Identification of Metabolites
#
# CMakeLists.txt
#
# Copyright (c) 2020, Fei Wang
# All rights reserved.
#
# This file is part of the cfm-id project.
# The contents are covered by the terms of the GNU Lesser General Public
# License, which is included in the file license.txt, found at the root
# of the cfm source tree.
#########################################################################*/
# We may need higher verison cmake
cmake_minimum_required(VERSION 3.7.0)
set(CMAKE_CXX_STANDARD 11)

# Set a default build type if none was specified
set(default_build_type "RelWithDebInfo")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    set(default_build_type "Debug")
endif()

project(cfm VERSION 4.0.6)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Version.h.in ${CMAKE_CURRENT_SOURCE_DIR}/cfm-code/Version.h)

set(CFM_OUTPUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE STRING "Where do you want the output binaries to install to?")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#It is recommended, but not strictly necessary, to build and run the tests
set(INCLUDE_TESTS OFF CACHE BOOL "Do you want to build the tests?")

#Building the training code requires extra dependencies, e.g. MPI
set(INCLUDE_TRAIN OFF CACHE BOOL "Do you want to build the training code?")

# BOOST
set(Boost_USE_STATIC_LIBS OFF CACHE BOOL "use static libraries from Boost")
include(FindBoost)

set(Boost_componets filesystem system serialization program_options)
# For some reason msvc will give LNK2019 if not include thread
if (MSVC)
    list(APPEND Boost_componets thread)
endif()

if (INCLUDE_TESTS)
    list(APPEND Boost_componets unit_test_framework)
endif()

# we need 1.62 or higher for unit testing
# on paper 1.59 + should work, seems 1.60 has building problem
find_package(Boost 1.62 REQUIRED COMPONENTS ${Boost_componets})

include_directories(${Boost_INCLUDE_DIRS})

#include cmake finds
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# RDKit
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR})
include(FindRDKit)

if (RDKIT_FOUND)
    include_directories(${RDKIT_INCLUDE_DIR} ${RDKIT_INCLUDE_EXT_DIR})
    link_directories(${RDKIT_LIBRARY_DIR})
else ()
    message(FATAL_ERROR "Could not find RDKit")
endif ()

# LPSolve
set(LPSOLVE_INCLUDE_DIR "" CACHE STRING "Where are the LPSolve headers?")
set(LPSOLVE_LIBRARY_DIR "" CACHE STRING "Where are the LPSolve libraries?")
message(STATUS "Search For lPsolve lib at ${LPSOLVE_LIBRARY_DIR} and headers at ${LPSOLVE_INCLUDE_DIR} ")


find_path(LPSOLVE_LIBRARY_HEADERS lp_lib.h PATHS
                        ${LPSOLVE_INCLUDE_DIR}
                        /usr/include/lpsolve
                        /usr/local/include/lpsolve
                        /opt/lpsolve
                        ~/lpsolve
                        /usr/include/lp_solve
                        /usr/local/include/lp_solve
                        /opt/lp_solve
                        ~/lp_solve)

if (NOT LPSOLVE_LIBRARY_HEADERS)
    MESSAGE(FATAL_ERROR "Cannot find LPSOLVE Headers!")
else()
    MESSAGE(STATUS "Found LPSOLVE Headers at ${LPSOLVE_LIBRARY_HEADERS}")
endif (NOT LPSOLVE_LIBRARY_HEADERS)


 find_library(LPSOLVE_LIBRARY lpsolve55 PATHS
                         ${LPSOLVE_LIBRARY_DIR}
                        /usr/include/lpsolve
                        /usr/local/include/lpsolve
                        /opt/lpsolve
                        ~/lpsolve
                        /usr/include/lp_solve
                        /usr/local/include/lp_solve
                        /opt/lp_solve
                        ~/lp_solve)

if (NOT LPSOLVE_LIBRARY)
    MESSAGE(FATAL_ERROR "Cannot find LPSOLVE libraries!")
else()
    MESSAGE(STATUS "Found LPSOLVE libraries file ${LPSOLVE_LIBRARY}")
endif (NOT LPSOLVE_LIBRARY)

include_directories(${LPSOLVE_LIBRARY_HEADERS})
link_directories(LPSOLVE_LIBRARY)
mark_as_advanced(LPSOLVE_LIBRARY)

find_package(Threads)

#  MPI
if (INCLUDE_TRAIN OR INCLUDE_TESTS)
    include(FindMPI)
    if (MPI_FOUND)
        include_directories(${MPI_INCLUDE_PATH})
    else ()
        message(FATAL_ERROR "Could not find MPI")
    endif ()
endif ()

#Common Library
add_subdirectory(cfm-code)

#Executable for generating fragmentation graphs
include_directories(${CMAKE_MODULE_PATH}/cfm-code)
add_subdirectory(fraggraph-gen)

#Executable for running tests
if (INCLUDE_TESTS)
    add_subdirectory(cfm-test)
    enable_testing()
    add_subdirectory(cfm-boost-test)
endif ()

if (INCLUDE_TRAIN)
    add_subdirectory(cfm-train)
endif ()

#Spectrum Predicter
add_subdirectory(cfm-predict)

#Identifier
add_subdirectory(cfm-id)

#Annotater
add_subdirectory(cfm-annotate)

#Compute Statistics (for computed spectra)
add_subdirectory(compute-stats)

add_subdirectory(cfm-id-precomputed)