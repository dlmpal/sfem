#pragma once

#include "../../kernel.h"

namespace sfem::kernel::elasticity
{

    static int const PLANE_STRESS = 0;
    static int const PLANE_STRAIN = 1;
    static int const SOLID = 2;

    class Elasticity : public Kernel
    {
    public:
        Elasticity(Float E, Float nu);
        virtual ~Elasticity() = 0;

    protected:
        virtual std::pair<int, int> ShapeDerivMatrixSize() = 0;
        virtual void ShapeDerivMatrix(Float B[]) = 0;
        virtual void Constitutive(Float D[]) = 0;
        void operator()(Float kloc[]) override;

        /* Young's Modulus */
        Float E;

        /* Poisson's ratio */
        Float nu;
    };
}