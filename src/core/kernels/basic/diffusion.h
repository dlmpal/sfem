#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{

    /// @brief Diffusion Kernel
    class Diffusion : public Kernel
    {
    public:
        Diffusion(Float c);
        KernelType GetType() const override;
        Float GetCoefficient() const;

    private:
        void Evaluate(std::vector<Float> &kloc) override;

        /// @brief Diffusion coefficient
        Float c;
    };
}