#==============================================================================
# Top level CMakeLists.txt for sfem
#==============================================================================
cmake_minimum_required(VERSION 3.16)
#==============================================================================
# Project name and version
project(SFEM VERSION "0.0.1")
#==============================================================================
# Include cmake configuration
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake)
#==============================================================================
# Build the library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
#==============================================================================
# Build applications
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/apps)