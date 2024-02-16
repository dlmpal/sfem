#==============================================================================
# Create the library
add_library(sfem SHARED)
#==============================================================================
# Subdirectories
set(SFEM_DIRS 
        common
        mesh
        geo
        io
        la
        fe
        core)

foreach(dir ${SFEM_DIRS})
        add_subdirectory(${dir})
endforeach()
#==============================================================================
target_include_directories(sfem PUBLIC ${SFEM_SOURCE_DIR})