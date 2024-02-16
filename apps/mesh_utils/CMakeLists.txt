#==============================================================================
include("../../cmake/config.cmake")
#=============================================================================
# gmshToSfem
add_executable(gmshToSfem ${CMAKE_CURRENT_SOURCE_DIR}/gmsh_to_sfem.cc)
target_include_directories(gmshToSfem PUBLIC "${SFEM_SOURCE_DIR}")
target_link_directories(gmshToSfem PUBLIC "${SFEM_BUILD_DIR}")
target_link_libraries(gmshToSfem sfem)
#==============================================================================
# sfemToVTK
add_executable(sfemToVTK ${CMAKE_CURRENT_SOURCE_DIR}/sfem_to_vtk.cc)
target_include_directories(sfemToVTK PUBLIC "${SFEM_SOURCE_DIR}")
target_link_directories(sfemToVTK PUBLIC "${SFEM_BUILD_DIR}")
target_link_libraries(sfemToVTK sfem)
#==============================================================================
# sfemPartMesh
add_executable(sfemPartMesh ${CMAKE_CURRENT_SOURCE_DIR}/sfem_part_mesh.cc)
target_include_directories(sfemPartMesh PUBLIC "${SFEM_SOURCE_DIR}")
target_link_directories(sfemPartMesh PUBLIC "${SFEM_BUILD_DIR}")
target_link_libraries(sfemPartMesh sfem)