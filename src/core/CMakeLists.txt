target_sources(sfem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/field.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/kernel.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/mono_field_assembler.cc)

add_subdirectory(kernels)
add_subdirectory(time_integrators)
add_subdirectory(solvers)