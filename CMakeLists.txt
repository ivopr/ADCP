cmake_minimum_required(VERSION 3.16)
project(ADCP VERSION 0.1 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(ADCP_NO_COMMON "Use -fno-common to catch symbol conflicts (recommended for CI)" OFF)

if(ADCP_NO_COMMON)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-common")
else()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fcommon")
endif()

option(ADCP_MPI    "Build with MPI parallel tempering support" OFF)
option(ADCP_OPENMP "Build with OpenMP support (required for cdlearn)" OFF)
option(ADCP_TOOLS  "Build auxiliary tools" ON)
option(ADCP_LEGACY "Build legacy metropolis_old binary" OFF)

if(ADCP_MPI)
  find_package(MPI REQUIRED)
endif()

if(ADCP_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

add_subdirectory(src)
if(ADCP_TOOLS)
  add_subdirectory(tools)
endif()

option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(tests)
endif()
