#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{
    class Source : public Kernel
    {
    public:
        Source(const Float c[]);

    private:
        virtual void ComputeValues(Float values[]);
        void operator()(Float kloc[]) override;

        /// @brief Source value for every variable
        const Float *c;
    };
}