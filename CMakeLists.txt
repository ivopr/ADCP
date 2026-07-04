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

  # ── Internal: run ctest without rebuild ──
  add_custom_target(run-ctest
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
            --output-on-failure -j4
    USES_TERMINAL
  )

  # ── check: rebuild + run all available tests ──
  set(_check_deps adcp cdlearn cm rama lipa mergie statistics dssp2cm oops)
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

  # ── check-phys / check-all: only if physical validation is opted in ──
  if(ADCP_VALIDATION)
    add_custom_target(check-phys
      COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
              -L validation --output-on-failure
      USES_TERMINAL
      COMMENT "Running ADCP physical validation tests"
    )

    add_custom_target(check-all
      COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR}
              --output-on-failure -j4
      DEPENDS ${_check_deps}
      USES_TERMINAL
      COMMENT "Running ALL ADCP tests (smoke + functional + validation)"
    )
  endif()
endif()
