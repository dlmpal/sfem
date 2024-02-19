#pragma once

#include "../../kernel.h"

namespace sfem::kernel::elasticity
{
    class Elasticity : public Kernel
    {
    public:
        Elasticity(Float E, Float nu);
        KernelType GetType() const override;
        virtual ~Elasticity() = 0;

    protected:
        virtual std::pair<int, int> ShapeDerivMatrixSize() = 0;
        virtual void ShapeDerivMatrix(Float B[]) = 0;
        virtual void Constitutive(Float D[]) = 0;
        void Evaluate(std::vector<Float> &kloc) override;

        /// @brief Young's Modulus
        Float E;

        /// @brief Poisson's ratio
        Float nu;
    };
}