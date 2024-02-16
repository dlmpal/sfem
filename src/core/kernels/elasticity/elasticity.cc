#include "elasticity.h"
#include "../../../common/mat_ops.h"

namespace sfem::kernel::elasticity
{
    //=============================================================================
    Elasticity::Elasticity(Float E, Float nu) : Kernel(KernelType::LHS)
    {
        this->E = E;
        this->nu = nu;
    }
    //=============================================================================
    Elasticity::~Elasticity()
    {
    }
    //=============================================================================
    void Elasticity::operator()(Float kloc[])
    {
        auto [n_rows, n_cols] = ShapeDerivMatrixSize();
        Float D[n_rows * n_rows] = {0.0};
        Float B[n_rows * n_cols] = {0.0};
        Float B_trans[n_cols * n_rows] = {0.0};
        Float D_B[n_rows * n_cols] = {0.0};

        ShapeDerivMatrix(B);
        mat_ops::mat_transpose(n_rows, n_cols, B, B_trans);
        Constitutive(D);
        mat_ops::mat_mult(n_rows, n_cols, n_rows, D, B, D_B);
        mat_ops::mat_mult(n_cols, n_cols, n_rows, B_trans, D_B, kloc);
    }

}