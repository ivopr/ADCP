target_include_directories(adcp_core PUBLIC ${CMAKE_SOURCE_DIR}/include)

add_executable(cdlearn cdlearn.cpp)
target_include_directories(cdlearn PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(cdlearn PRIVATE adcp_core)
if(ADCP_OPENMP)
  target_link_libraries(cdlearn PRIVATE OpenMP::OpenMP_CXX)
endif()

add_executable(cm cm.cpp)
target_include_directories(cm PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(cm PRIVATE adcp_core)

add_executable(rama ramachandran.cpp)
target_include_directories(rama PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(rama PRIVATE adcp_core)

add_executable(lipa pauling.cpp)
target_include_directories(lipa PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(lipa PRIVATE adcp_core)

add_executable(bfactor bfactor.cpp)
target_include_directories(bfactor PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(bfactor PRIVATE adcp_core)

add_executable(mergie mergie.cpp)
target_link_libraries(mergie PRIVATE m)

add_executable(statistics statistics.cpp)
target_link_libraries(statistics PRIVATE m)

add_executable(dssp2cm dssp2cm.cpp)
target_link_libraries(dssp2cm PRIVATE m)

add_executable(oops oops.cpp)
target_link_libraries(oops PRIVATE m)

install(TARGETS cdlearn cm rama lipa bfactor mergie statistics dssp2cm oops
  RUNTIME DESTINATION bin
)
