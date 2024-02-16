cmake_minimum_required(VERSION 3.16)
project(sfem_laplace_solver)

include(../../cmake/config.cmake)

# laplaceSolver
add_executable(laplaceSolver ../main.cc)
target_include_directories(laplaceSolver PUBLIC ${SFEM_SOURCE_DIR})
target_link_directories(laplaceSolver PUBLIC ${SFEM_BUILD_DIR})
target_link_libraries(laplaceSolver PUBLIC sfem)