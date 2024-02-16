#pragma once

#include "../../kernel.h"

namespace sfem::kernel::elasticity
{
    class Gravity : public Kernel
    {
    public:
        Gravity(Float rho, int direction);

    private:
        void operator()(Float kloc[]) override;

        /// @brief Gravitational acceleration
        Float g = 9.81;

        /// @brief Material density
        Float rho;

        /// @brief Direction in which gravity acts
        int direction;
    };
}
