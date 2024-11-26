cmake_minimum_required(VERSION 3.16)

option(enable_OPENMP "Enable OpenMP" ON)
option(enable_GLOG "Enable GLOG" OFF)

#-----compiler options
if(NOT DEFINED compiler)
  set(base_compile_options "-Wall -Wextra -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -fopenmp")
   endif()
endif()

if(compiler STREQUAL gcc)
  set(CMAKE_C_COMPILER "mpicc")
  set(CMAKE_CXX_COMPILER "mpicxx")
  set(base_compile_options "-Wall -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -fopenmp")
  endif()
elseif(DEFINED compiler)
  message(FATAL ERROR "invalid compiler")
endif()

if(compiler STREQUAL intel)
  set(CMAKE_C_COMPILER "mpiicx")
  set(CMAKE_CXX_COMPILER "mpiicpx")
  set(base_compile_options "-Wall -Wextra -O3")
  if(enable_OPENMP)
    set(base_compile_options ${base_compile_options} " -qopenmp")
  endif()
elseif(DEFINED compiler)
  message(FATAL ERROR "invalid compiler")
endif()
string(REPLACE " " ";" compile_options ${base_compile_options})
#-----------------------

project(topologyOptim C CXX)
set(CMAKE_CXX_STANDARD 17)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

#change default install prefix
set(INSTALL_DIR "/root/c++/AdvectionDiffusionSolver")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${INSTALL_DIR} CACHE PATH "..." FORCE)
endif()

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)

#######
# Display options
#######
message(" ")
message("${PROJECT_NAME}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "COMPILE_OPTIONS       : "      ${compile_options})
message(" ")

#######
# find Libraries to depend on
###


# Google-glog
if(enable_GLOG)
  set(GLOG_DIR ${GLOG_DIR} CACHE PATH "...")
  find_package(GLOG REQUIRED)
  include_directories(${GLOG_INCLUDE_DIR})
  message( STATUS "GLOG include dir  : "      ${GLOG_INCLUDE_DIR})
  MESSAGE("find GLOG")
endif()

# Eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN_INCLUDE_DIR})

#######
# Sub Directories
#######

# add_subdirectory(include)
add_subdirectory(src)