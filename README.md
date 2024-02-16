# 𝙎𝙁𝙀𝙈

SFEM is a finite element toolbox written in C++. It supports both sequentual and parallel computation,
by utilizing the MPI protocol and PETSc. Parallel execution also reqires METIS for mesh partitioning.
SLEPc is also an optional requirement for eigenvalue computation of discretized operators.

## Requirements
* C++ 17 and above
* PETSc (linear algebra)
* MPI (optional, parallel execution)
* METIS (optional, mesh partioning)
* SLEPc (optional, eigenvalue computation)

## Installation
Start by cloning the repository to your machine. Open the config.cmake file, found under cmake/.
There edit the path your PETSc and MPI installations. The same can be done for any other optional dependency 
(remember to also turn on the corresponding flag, e.g SFEM_USE_SLEPC). Finally, while in the SFEM home directory
create a build folder. Once inside it, run cmake and then make.

## Example Programs
Example use of SFEM can be found under the tests/ folder.

## Non-native Mesh Formats
Currently the only non-native formats supported by SFEM are Gmsh for importing,
and VTK for exporting meshes. Programs for doing so are found under apps/mesh_utils.