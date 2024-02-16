# Basic kernels
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/basic/diffusion.cc
${CMAKE_CURRENT_SOURCE_DIR}/basic/mass.cc
${CMAKE_CURRENT_SOURCE_DIR}/basic/source.cc)

# Elasiticty kernels
target_sources(sfem PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/elasticity/elasticity.cc
${CMAKE_CURRENT_SOURCE_DIR}/elasticity/elasticity_2d.cc
${CMAKE_CURRENT_SOURCE_DIR}/elasticity/elasticity_3d.cc
${CMAKE_CURRENT_SOURCE_DIR}/elasticity/gravity.cc
${CMAKE_CURRENT_SOURCE_DIR}/elasticity/pressure.cc
)