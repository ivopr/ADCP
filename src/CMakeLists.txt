add_library(adcp_core STATIC
  aadict.c
  canonicalAA.c
  checkpoint_io.c
  energy.c
  error.c
  flex.c
  metropolis.c
  nested.c
  params.c
  peptide.c
  probe.c
  random16.c
  rotation.c
  vdw.c
  vector.c
)
target_include_directories(adcp_core PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)
target_link_libraries(adcp_core PUBLIC m)

add_executable(adcp main.c)
target_link_libraries(adcp PRIVATE adcp_core)

if(ADCP_MPI)
  add_executable(adcp_mpi main.c)
  target_compile_definitions(adcp_mpi PRIVATE PARALLEL)
  target_include_directories(adcp_mpi PRIVATE ${MPI_C_INCLUDE_DIRS})
  target_link_libraries(adcp_mpi PRIVATE adcp_core ${MPI_C_LIBRARIES})
endif()

install(TARGETS adcp adcp_core
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
)
