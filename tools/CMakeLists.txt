target_include_directories(adcp_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_executable(cdlearn cdlearn.c)
target_include_directories(cdlearn PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(cdlearn PRIVATE adcp_core)
if(ADCP_OPENMP)
  target_link_libraries(cdlearn PRIVATE OpenMP::OpenMP_C)
endif()

add_executable(cm cm.c)
target_include_directories(cm PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(cm PRIVATE adcp_core)

add_executable(rama ramachandran.c)
target_include_directories(rama PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(rama PRIVATE adcp_core)

add_executable(lipa pauling.c)
target_include_directories(lipa PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lipa PRIVATE adcp_core)

add_executable(bfactor bfactor.c)
target_include_directories(bfactor PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(bfactor PRIVATE adcp_core)

add_executable(mergie mergie.c)
target_link_libraries(mergie PRIVATE m)

add_executable(statistics statistics.c)
target_link_libraries(statistics PRIVATE m)

add_executable(dssp2cm dssp2cm.c)
target_link_libraries(dssp2cm PRIVATE m)

add_executable(oops oops.c)
target_link_libraries(oops PRIVATE m)

if(ADCP_LEGACY)
  add_executable(metropolis_old metropolis_old.c)
  target_include_directories(metropolis_old PRIVATE ${CMAKE_SOURCE_DIR}/include)
  target_link_libraries(metropolis_old PRIVATE adcp_core)
endif()

install(TARGETS cdlearn cm rama lipa bfactor mergie statistics dssp2cm oops
  RUNTIME DESTINATION bin
)
