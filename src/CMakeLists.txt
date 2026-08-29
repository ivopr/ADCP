set(ADCP_CORE_SOURCES
  aadict.cpp
  canonicalAA.cpp
  checkpoint_io.cpp
  energy.cpp
  error.cpp
  flex.cpp
  metropolis.cpp
  nested.cpp
  params.cpp
  peptide.cpp
  probe.cpp
  random16.cpp
  rotation.cpp
  vdw.cpp
  vector.cpp
)

add_library(adcp_core STATIC ${ADCP_CORE_SOURCES})
target_include_directories(adcp_core PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(adcp_core PUBLIC m)

add_executable(adcp main.cpp)
target_link_libraries(adcp PRIVATE adcp_core)

if(ADCP_MPI)
  # adcp_core above is compiled without PARALLEL, so nested.cpp/checkpoint_io.cpp/
  # flex.cpp/probe.cpp never see the MPI code paths. A second static lib, compiled
  # from the same sources with PARALLEL defined, is what actually builds them.
  add_library(adcp_core_mpi STATIC ${ADCP_CORE_SOURCES})
  target_include_directories(adcp_core_mpi PUBLIC
    ${CMAKE_SOURCE_DIR}/include
  )
  target_compile_definitions(adcp_core_mpi PUBLIC PARALLEL)
  target_link_libraries(adcp_core_mpi PUBLIC m MPI::MPI_CXX)

  add_executable(adcp_mpi main.cpp)
  target_compile_definitions(adcp_mpi PRIVATE PARALLEL)
  target_link_libraries(adcp_mpi PRIVATE adcp_core_mpi MPI::MPI_CXX)
endif()

install(TARGETS adcp adcp_core
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
)
