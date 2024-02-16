#include "thermal_stress.h"

namespace sfem::kernel::elasticity
{
    ThermalStress::ThermalStress(Float alpha, Float dT, Elasticity *kernel) : Kernel(ADD_RHS)
    {
        this->alpha = alpha;

        this->dT = dT;

        this->kernel = kernel;
    }

    void ThermalStress::operator()(Float kloc[])
    {

        auto size = kernel->GetShapeDerivMatrixSize();

        auto n_rows = size.first;

        auto n_cols = size.second;

        Float thermal_strain[n_rows] = {0.0};

        Float coeff = alpha * dT;

        if (kernel->analysis_type == PLANE_STRAIN)

            coeff *= (1 + kernel->nu);

        for (auto i = 0; i < fe.dim; i++)

            thermal_strain[i] = coeff;

        Float D[n_rows * n_rows] = {0.0};

        Float B[n_rows * n_cols] = {0.0};

        Float B_trans[n_cols * n_rows] = {0.0};

        Float B_trans_D[n_cols * n_rows] = {0.0};

        kernel->ShapeDerivMatrix(fe, B);

        mat_ops::mat_transpose(n_rows, n_cols, B, B_trans);

        kernel->GetConstitutive(D);

        mat_ops::mat_mult(n_cols, n_rows, n_rows, B_trans, D, B_trans_D);

        for (auto i = 0; i < n_cols; i++)

            for (auto j = 0; j < n_rows; j++)

                kloc[i] += B_trans_D[i * n_rows + j] * thermal_strain[j];
    }
}