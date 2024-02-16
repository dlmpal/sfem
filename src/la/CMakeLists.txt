target_sources(sfem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sparsity_pattern.cc)

if(SFEM_USE_PETSC)
    target_sources(sfem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/petsc.cc)
endif()

if(SFEM_USE_SLEPC)
    target_sources(sfem PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/slepc.cc)
endif()