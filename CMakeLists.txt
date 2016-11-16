
add_library(KMCode "compressed_mesh" )
target_link_libraries(mmpl_adapter ddl_parmetis ${mmpd_ADAPTER_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES COMPILE_DEFINITIONS "${MODFEM_MPI_CODE_DEFINE}")
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

enable_testing()

add_executable(MOD_FEM_test_mmpl_adapter "mmps_adapter_test.cpp")
target_link_libraries(MOD_FEM_test_mmpl_adapter utl_log utl_io mml_t4_prism pcl_mpi_safe mmpl_adapter  ${Boost_LIBRARIES})

add_test(NAME MOD_FEM_test_mmpl_adapter
COMMAND mpirun -np 2 ${CMAKE_BINARY_DIR}/MOD_FEM_test_mmpl_adapter
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")

