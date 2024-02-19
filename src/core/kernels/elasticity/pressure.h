#pragma once

#include "../../kernel.h"
#include "../../field.h"

namespace sfem::kernel::elasticity
{
    class Pressure : public Kernel
    {
    public:
        Pressure(field::Field *p);
        KernelType GetType() const override;

    private:
        void Evaluate(std::vector<Float> &kloc);

        /// @brief Pressure field
        field::Field *p;
    };
}