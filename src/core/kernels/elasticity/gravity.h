#pragma once

#include "../../kernel.h"

namespace sfem::kernel::elasticity
{
    class Gravity : public Kernel
    {
    public:
        Gravity(Float rho, int direction);
        KernelType GetType() const override;

    private:
        void Evaluate(std::vector<Float> &kloc) override;

        /// @brief Gravitational acceleration
        Float g = 9.81;

        /// @brief Material density
        Float rho;

        /// @brief Direction in which gravity acts
        int direction;
    };
}
