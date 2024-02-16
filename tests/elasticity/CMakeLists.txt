cmake_minimum_required(VERSION 3.16)
project(sfem_elasticity_solver)

include(../../cmake/config.cmake)

# elasticitySolver
add_executable(elasticitySolver ../main.cc)
target_include_directories(elasticitySolver PUBLIC ${SFEM_SOURCE_DIR})
target_link_directories(elasticitySolver PUBLIC ${SFEM_BUILD_DIR})
target_link_libraries(elasticitySolver sfem)