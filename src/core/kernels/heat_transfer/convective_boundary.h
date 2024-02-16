#pragma once

#include "../basic/diffusion.h"

namespace sfem::kernel::heat_transfer
{
    class ConvectiveBoundary : public Kernel
    {
    public:
        ConvectiveBoundary(Float h, Float T_inf, basic::Diffusion *kernel);

    private:
        void operator()(Float kloc[]) override;

        /* Convection coefficient */
        Float h;

        /* Free-stream temperature */
        Float T_inf;

        /* Heat conduction kernel */
        basic::Diffusion *kernel;
    };
}