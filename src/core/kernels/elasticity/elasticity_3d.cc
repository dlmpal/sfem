#include "elasticity_3d.h"

namespace sfem::kernel::elasticity
{
    //=============================================================================
    Elasticity3D::Elasticity3D(Float E, Float nu) : Elasticity(E, nu)
    {
    }
    //=============================================================================
    std::pair<int, int> Elasticity3D::ShapeDerivMatrixSize()
    {
        int n_rows = 6;
        int n_cols = 3 * fe->GetNumNodes();
        return std::make_pair(n_rows, n_cols);
    }
    //=============================================================================
    void Elasticity3D::ShapeDerivMatrix(Float B[])
    {
        auto [n_rows, n_cols] = ShapeDerivMatrixSize();
        auto dNdX = fe->dNdX();
        for (auto i = 0; i < fe->GetNumNodes(); i++)
        {
            /* exx */
            B[0 * n_cols + i * 3 + 0] = dNdX[i][0];

            /* eyy */
            B[1 * n_cols + i * 3 + 1] = dNdX[i][1];

            /* ezz */
            B[2 * n_cols + i * 3 + 2] = dNdX[i][2];

            /* exy */
            B[3 * n_cols + i * 3 + 0] = dNdX[i][1];
            B[3 * n_cols + i * 3 + 1] = dNdX[i][0];

            /* eyz */
            B[4 * n_cols + i * 3 + 1] = dNdX[i][2];
            B[4 * n_cols + i * 3 + 2] = dNdX[i][1];

            /* exz */
            B[5 * n_cols + i * 3 + 0] = dNdX[i][2];
            B[5 * n_cols + i * 3 + 2] = dNdX[i][0];
        }
    }
    //=============================================================================
    void Elasticity3D::Constitutive(Float D[])
    {
        Float c1 = E / ((1 + nu) * (1 - 2 * nu));
        Float c2 = (1 - 2 * nu) / 2;

        D[0] = (1 - nu) * c1;

        D[1] = nu * c1;

        D[2] = nu * c1;

        D[6] = nu * c1;

        D[7] = (1 - nu) * c1;

        D[8] = nu * c1;

        D[12] = nu * c1;

        D[13] = nu * c1;

        D[14] = (1 - nu) * c1;

        D[21] = c1 * c2;

        D[28] = c1 * c2;

        D[35] = c1 * c2;
    }
}