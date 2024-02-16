#==============================================================================
# C++ standard
set(CMAKE_CXX_STANDARD 17)
#==============================================================================
# Enforce standard
set(CMAKE_CXX_STANDARD_REQUIRED on)
#==============================================================================
# Build type
## Relase (-O3)
if (CMAKE_BUILD_TYPE MATCHES RELEASE)
        set(CMAKE_CXX_FLAGS "-O3")
endif (CMAKE_BUILD_TYPE MATCHES RELEASE)

## Debug
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
        set(CMAKE_CXX_FLAGS "-g -Wall")
endif (CMAKE_BUILD_TYPE MATCHES DEBUG)
#==============================================================================
# SFEM home, source and build directory
set(SFEM_DIR ${CMAKE_CURRENT_LIST_DIR}/..)
set(SFEM_SOURCE_DIR ${SFEM_DIR}/src)
set(SFEM_BUILD_DIR ${SFEM_DIR}/build/src)
#==============================================================================
# Compile definitions used by SFEM
set(SFEM_COMPILE_DEFINITIONS "")
#==============================================================================
# Third-party dependencies
# The following lists contain the include and link directories
# and the libraries used by SFEM and its dependencies.
# They are added to all targets including this file.
set(SFEM_INCLUDE_DEPS "")
set(SFEM_LINK_DEPS "")
set(SFEM_LIB_DEPS "")

# PETSc
set(SFEM_USE_PETSC on)
if(SFEM_USE_PETSC)
    set(PETSC_DIR /usr/lib/petsc)
    list(APPEND SFEM_INCLUDE_DEPS ${PETSC_DIR}/include)
    list(APPEND SFEM_LINK_DEPS ${PETSC_DIR}/lib)
    list(APPEND SFEM_LIB_DEPS petsc_real)
    list(APPEND SFEM_COMPILE_DEFINITIONS SFEM_USE_PETSC)
else()
    message(FATAL_ERROR "SFEM currently only supports PETSc as the linear algebra backend!")
endif()

## SLEPc
set(SFEM_USE_SLEPC off)
if(SFEM_USE_SLEPC)
if(NOT SFEM_USE_PETSC)
    message(FATAL_ERROR "SLEPc requires PETSc backend!")
endif()
    set(SLEPC_DIR /usr/lib/slepc)
    list(APPEND SFEM_INCLUDE_DEPS ${SLEPC_DIR}/include)
    list(APPEND SFEM_LINK_DEPS ${SLEPC_DIR}/lib)
    list(APPEND SFEM_LIB_DEPS slepc_real)
    list(APPEND SFEM_COMPILE_DEFINITIONS SFEM_USE_SLEPC)
endif()

## METIS
set(SFEM_USE_METIS off)
if(SFEM_USE_METIS)
    set(METIS_DIR /home/user/local)
    list(APPEND SFEM_INCLUDE_DEPS ${METIS_DIR}/include)
    list(APPEND SFEM_LINK_DEPS ${METIS_DIR}/lib)
    list(APPEND SFEM_LIB_DEPS metis)
    list(APPEND SFEM_COMPILE_DEFINITIONS SFEM_USE_METIS)
endif()

## MPI
set(SFEM_USE_MPI on)
if(SFEM_USE_MPI)
    set(MPI_DIR /usr/lib/x86_64-linux-gnu/openmpi)
    list(APPEND SFEM_INCLUDE_DEPS ${MPI_DIR}/include)
    list(APPEND SFEM_LINK_DEPS ${MPI_DIR}/lib)
    list(APPEND SFEM_COMPILE_DEFINITIONS SFEM_USE_MPI)
endif()

# Add directories and libs to targets
include_directories(${SFEM_INCLUDE_DEPS})
link_directories(${SFEM_LINK_DEPS})
link_libraries(${SFEM_LIB_DEPS})
#==============================================================================
# Add compile definitions to targets
add_compile_definitions(${SFEM_COMPILE_DEFINITIONS})
#==============================================================================
# Select compiler
if(SFEM_USE_MPI)
    set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)     
else()
    set(CMAKE_CXX_COMPILER /usr/bin/g++)
endif()
#==============================================================================
