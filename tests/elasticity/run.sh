#!/bin/bash

mkdir -p build
rm -rf build/*

(cd ../../build; make)
(cd build; cmake -S .. -B .; make)

mkdir -p fields
rm fields/*

# Modal solver
# mpiexec -np 4 build/elasticitySolver -st_type sinvert
# (cd fields; ../../../applications/mesh_utilities/build/sfemToVTK ../Mesh3D . 1 u_mode_0 u_mode_1 u_mode_2)

# Static solver
mpiexec -np 4 build/elasticitySolver -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package_mumps
(cd fields; ../../../build/apps/mesh_utils/sfemToVTK ../Mesh3D 1 U)
paraview fields/sfem_0.vtk