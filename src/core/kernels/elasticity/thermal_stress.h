#pragma once

#include "elasticity.h"

namespace sfem::kernel::elasticity
{
    struct ThermalStress : Kernel
    {

        /* Thermal expansion coefficient */
        Float alpha;

        /* Temperature difference */
        Float dT;

        /* Corresponding elasticity kernel  */
        Elasticity *kernel;

        ThermalStress(Float alpha, Float dT, Elasticity *kernel);

        void operator()(Float kloc[]) override;
    };
}