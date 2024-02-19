#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{
    /// @brief Mass Kernel
    class Mass : public Kernel
    {
    public:
        Mass(Float c);
        KernelType GetType() const override;
        Float GetCoefficient() const;

    private:
        void Evaluate(std::vector<Float> &kloc) override;

        /// @brief Mass matrix coefficient
        Float c;
    };
}
