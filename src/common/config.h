#pragma once

// C++ std-lib headers
#include "utility"
#include "string"
#include "vector"
#include "array"
#include "map"
#include "unordered_map"
#include "memory"
#include "algorithm"
#include "functional"
#include "cmath"
#include "iostream"
#include "fstream"

// Data type used for floating-point arithmetic
#define Float double

// Root process rank
#define SFEM_ROOT 0

// Max number of nodes per cell
#define SFEM_MAX_CELL_NODES 30

// Max number of variables per field
#define SFEM_MAX_N_VARS_FIELD 6

// Max number of quadrature points
#define SFEM_MAX_QUAD_POINTS 10

// MPI
#ifdef SFEM_USE_MPI
#include "mpi.h"

// Global MPI communicator
#define SFEM_COMM_WORLD MPI_COMM_WORLD

// Float datatype used for MPI operations
#define SFEM_MPI_FLOAT MPI_DOUBLE
#endif // FEMXX_USE_MPI

namespace sfem
{

}