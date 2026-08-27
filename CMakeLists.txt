cmake_minimum_required(VERSION 3.16)
project(ADCP VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ADCP_LEGACY_COMMON "Allow tentative global merging (-fcommon); legacy escape hatch" OFF)

if(ADCP_LEGACY_COMMON)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcommon")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-common")
endif()

option(ADCP_MPI    "Build with MPI parallel tempering support" OFF)
option(ADCP_OPENMP "Build with OpenMP support (required for cdlearn)" OFF)
option(ADCP_TOOLS  "Build auxiliary tools" ON)

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

  # ── Internal: run ctest without rebuild ──
  add_custom_target(run-ctest
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
            --output-on-failure -j4
    USES_TERMINAL
  )

  # ── check: rebuild + run all available tests ──
  set(_check_deps adcp cdlearn cm rama lipa bfactor mergie statistics dssp2cm oops)
  if(NOT ADCP_TOOLS)
    set(_check_deps adcp)
  endif()

  add_custom_target(check
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
            --output-on-failure -j4
    DEPENDS ${_check_deps}
    USES_TERMINAL
    COMMENT "Running all ADCP tests (smoke + functional)"
  )

endif()
