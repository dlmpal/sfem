#pragma once

/// @brief Standard FE Kernels
namespace sfem::kernel
{

}

// Basic kernels
#include "basic/diffusion.h"
#include "basic/mass.h"
#include "basic/source.h"

// Elasticity kernels
#include "elasticity/sfem_elasticity.h"

// Heat transfer kernels
#include "heat_transfer/sfem_heat_transfer.h"