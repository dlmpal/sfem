#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{
    /// @brief Mass Kernel
    class Mass : public Kernel
    {
    public:
        Mass(Float c);
        Float GetCoefficient() const;

    private:
        void operator()(Float kloc[]) override;

        /// @brief Mass matrix coefficient
        Float c;
    };
}
