#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{

    /// @brief Diffusion Kernel
    class Diffusion : public Kernel
    {
    public:
        Diffusion(Float D);
        Float GetCoefficient() const;

    private:
        void operator()(Float kloc[]) override;

        /// @brief Diffusion coefficient
        Float c;
    };
}