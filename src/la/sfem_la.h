#pragma once

/// @brief Linear Algebra
namespace sfem::la
{

}

#include "sparsity_pattern.h"

#ifdef SFEM_USE_PETSC
#include "petsc.h"
#endif

#ifdef SFEM_USE_SLEPC
#include "slepc.h"
#endif