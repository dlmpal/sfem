#pragma once

#include "../../kernel.h"

namespace sfem::kernel::basic
{
    class Source : public Kernel
    {
    public:
        Source(const Float c[]);
        KernelType GetType() const override;

    private:
        virtual void ComputeValues(Float values[]);
        void Evaluate(std::vector<Float> &kloc) override;

        /// @brief Source value for every variable
        const Float *c;
    };
}