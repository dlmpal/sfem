#include "elasticity_2d.h"

namespace sfem::kernel::elasticity
{
    //=============================================================================
    Elasticity2D::Elasticity2D(Float E, Float nu) : Elasticity(E, nu)
    {
    }
    //=============================================================================
    std::pair<int, int> Elasticity2D::ShapeDerivMatrixSize()
    {
        int n_rows = 3;
        int n_cols = 2 * fe->GetNumNodes();
        return std::make_pair(n_rows, n_cols);
    }
    //=============================================================================
    void Elasticity2D::ShapeDerivMatrix(Float B[])
    {
        auto [n_rows, n_cols] = ShapeDerivMatrixSize();
        auto dNdX = fe->dNdX();
        for (auto i = 0; i < fe->GetNumNodes(); i++)
        {
            /* exx */
            B[0 * n_cols + i * 2 + 0] = dNdX[i][0];

            /* eyy */
            B[1 * n_cols + i * 2 + 1] = dNdX[i][1];

            /* exy */
            B[2 * n_cols + i * 2 + 0] = dNdX[i][1];
            B[2 * n_cols + i * 2 + 1] = dNdX[i][0];
        }
    }
    //=============================================================================
    PlaneStress::PlaneStress(Float E, Float nu) : Elasticity2D(E, nu)
    {
    }
    //=============================================================================
    void PlaneStress::Constitutive(Float D[])
    {
        Float coeff = E / (1 - nu * nu);

        D[0] = coeff * 1.0;
        D[1] = coeff * nu;
        D[2] = 0.0;

        D[3] = coeff * nu;
        D[4] = coeff * 1.0;
        D[5] = 0.0;

        D[6] = 0.0;
        D[7] = 0.0;
        D[8] = coeff * (1 - nu) * 0.5;
    }
    //=============================================================================
    PlaneStrain::PlaneStrain(Float E, Float nu) : Elasticity2D(E, nu)
    {
    }
    //=============================================================================
    void PlaneStrain::Constitutive(Float D[])
    {
        Float coeff = E / (1 - 2 * nu) / (1 + nu);

        D[0] = coeff * (1.0 - nu);
        D[1] = coeff * nu;
        D[2] = 0.0;

        D[3] = coeff * nu;
        D[4] = coeff * (1.0 - nu);
        D[5] = 0.0;

        D[6] = 0.0;
        D[7] = 0.0;
        D[8] = coeff * (1 - 2 * nu) * 0.5;
    }

}