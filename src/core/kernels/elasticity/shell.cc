#include "shell.h"

namespace sfem::kernel::elasticity
{

    ShellElement::ShellElement(const fe::basis::Basis *basis, const geo::Shape *shape, Float thick, const Float xpts[])
    {
        this->basis = basis;
        this->shape = shape;
        this->thick = thick;
        this->xpts = xpts;
        this->GetQuadrature();
        this->GetDirctors();
    }

    void ShellElement::GetQuadrature()
    {

        /*
         * The quadrature of the basis is augmented
         * with 2 points along the thickness
         */
        this->n_qpts = 2 * basis->n_qpts;
        for (auto i = 0; i < basis->n_qpts; i++)
        {
            for (auto j = 0; j < 2; j++)
            {
                qpts[i][j] = basis->qpts[i][j];
                qpts[i + basis->n_qpts][j] = basis->qpts[i][j];
            }
            qpts[i][2] = -0.577350269189626;
            qpts[i + basis->n_qpts][2] = 0.577350269189626;

            qWts[i] = basis->qwts[i];
            qWts[i + basis->n_qpts] = basis->qwts[i];
        }
    }

    void ShellElement::GetDirctors()
    {
        /*
         * v3 is the mid-surface normal,
         * therefore it is computed as the outward normal
         * of the triangle or quadrangle
         */
        v3 = shape->FaceNormal(-1);

        /*
         * v1 is perpendicular to v3.
         * min_x is the smallest component of v3
         * I is a unit vector in the direction of min_x
         * v1 is the cross product of I and v3 (v1 = I x v3)
         */
        Float min_x = v3.x[0];
        int min_i = 0;
        for (auto i = 1; i < 3; i++)
        {
            if (abs(v3.x[i]) < min_x)
            {
                min_x = abs(v3.x[i]);
                min_i = i;
            }
        }
        geo::Vec3 I(0, 0, 0);
        I.x[min_i] = 1.0;
        v1 = I.CrossProduct(v3);

        /*
         * v2 is perpendicular to both v1 and v3,
         * therefore v2 = v3 x v1
         */
        v2 = v3.CrossProduct(v1);

        G1 = v1 * (0.5 * thick);
        G2 = v2 * (-0.5 * thick);
    }

    void ShellElement::PhysicalToNaturalJac()
    {
        /*
         * x = N_k * (x_k + 0.5 * thick * v3 * zeta)
         */
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < basis->n_nodes; k++)
                {
                    dXdxi[i][j] += dNdxi[k][j] * xpts[k * 3 + i] + 0.5 * thick * v3.x[i] * pt[2];
                }
            }

            for (int k = 0; k < basis->n_nodes; k++)
            {
                dXdxi[i][2] += N[k] * 0.5 * thick * v3.x[i];
            }
        }
    }

    void ShellElement::NaturalToPhysicalJac()
    {
        J = mat_ops::inv(3, dXdxi, dxidX);
        if (J < 0)
        {
            std::string message = "Non-positive jacobian\n";
            Logger::GetInstance().LogMessage(message, Logger::WARN);
        }
    }

    void ShellElement::ComputeShapeGrad()
    {
        for (int i = 0; i < basis->n_nodes; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    dNdX[i][j] += dxidX[k][j] * dNdxi[i][k];
    }

    void ShellElement::ComputeTransform(const Float pt[])
    {
        for (auto i = 0; i < 3; i++)
        {
            for (auto j = 0; j < 3; j++)
            {
                dXdxi[i][j] = 0.0;
                dxidX[i][j] = 0.0;
            }
        }

        for (auto i = 0; i < basis->n_nodes; i++)
        {
            N[i] = 0.0;
            for (auto j = 0; j < 3; j++)
            {
                dNdxi[i][j] = 0.0;
                dNdX[i][j] = 0.0;
            }
        }

        J = -1;

        this->pt = pt;
        basis->ComputeShape(this->pt, N);
        basis->ComputeShapeAndGrad(this->pt, dNdxi);
        PhysicalToNaturalJac();
        NaturalToPhysicalJac();
        ComputeShapeGrad();
    }

    Shell::Shell(Float E, Float nu, Float thick) : Elasticity(E, nu)
    {
        this->thick = thick;
    }

    std::pair<int, int> Shell::ShapeDerivMatrixSize()
    {
        int n_rows = 5;
        int n_cols = 5 * fe->n_nodes;
        return std::make_pair(n_rows, n_cols);
    }

    void Shell::ShapeDerivMatrix(Float B[])
    {

        int n_rows = 5;
        int n_cols = 5 * shell->basis->n_nodes;
        mat_ops::mat_set(n_rows * n_cols, 0.0, B);
        Float dMidx, dMidy, dMidz;
        for (auto i = 0; i < shell->basis->n_nodes; i++)
        {
            dMidx = shell->dNdX[i][0] * shell->pt[2] + shell->N[i] * shell->dxidX[2][0];
            dMidy = shell->dNdX[i][1] * shell->pt[2] + shell->N[i] * shell->dxidX[2][1];
            dMidz = shell->dNdX[i][2] * shell->pt[2] + shell->N[i] * shell->dxidX[2][2];

            /* exx */
            B[0 * n_cols + i * 5 + 0] = shell->dNdX[i][0];
            B[0 * n_cols + i * 5 + 3] = dMidx * shell->G1.x[0];
            B[0 * n_cols + i * 5 + 4] = dMidx * shell->G2.x[0];

            /* eyy */
            B[1 * n_cols + i * 5 + 1] = shell->dNdX[i][1];
            B[1 * n_cols + i * 5 + 3] = dMidy * shell->G1.x[1];
            B[1 * n_cols + i * 5 + 4] = dMidy * shell->G2.x[1];

            /* exy */
            B[2 * n_cols + i * 5 + 0] = shell->dNdX[i][1];
            B[2 * n_cols + i * 5 + 1] = shell->dNdX[i][0];
            B[2 * n_cols + i * 5 + 3] = dMidx * shell->G1.x[1] + dMidy * shell->G1.x[0];
            B[2 * n_cols + i * 5 + 4] = dMidx * shell->G2.x[1] + dMidy * shell->G2.x[0];

            /* eyz */
            B[3 * n_cols + i * 5 + 1] = shell->dNdX[i][2];
            B[3 * n_cols + i * 5 + 2] = shell->dNdX[i][1];
            B[3 * n_cols + i * 5 + 3] = dMidy * shell->G1.x[2] + dMidz * shell->G1.x[1];
            B[3 * n_cols + i * 5 + 4] = dMidy * shell->G2.x[2] + dMidz * shell->G2.x[1];

            /* exz */
            B[4 * n_cols + i * 5 + 0] = shell->dNdX[i][2];
            B[4 * n_cols + i * 5 + 2] = shell->dNdX[i][0];
            B[4 * n_cols + i * 5 + 3] = dMidx * shell->G1.x[2] + dMidz * shell->G1.x[0];
            B[4 * n_cols + i * 5 + 4] = dMidx * shell->G2.x[2] + dMidz * shell->G2.x[0];
        }
        // mat_ops::mat_print(n_rows, n_cols, B);
        // exit(-1);
    }

    void Shell::Constitutive(Float D[])
    {
        Float k = 5.0 / 6.0;
        Float coeff = E / (1 - nu * nu);

        D[0] = coeff;
        D[1] = nu * coeff;
        D[2] = 0;
        D[3] = 0;
        D[4] = 0;

        D[5] = nu * coeff;
        D[6] = coeff;
        D[7] = 0;
        D[8] = 0;
        D[9] = 0;

        D[10] = 0;
        D[11] = 0;
        D[12] = 0.5 * (1 - nu) * coeff;
        D[13] = 0;
        D[14] = 0;

        D[15] = 0;
        D[16] = 0;
        D[17] = 0;
        D[18] = 0.5 * k * (1 - nu) * coeff;
        D[19] = 0;

        D[20] = 0;
        D[21] = 0;
        D[22] = 0;
        D[23] = 0;
        D[24] = 0.5 * k * (1 - nu) * coeff;
    }

    void Shell::Integrate(fe::FiniteElement *fe, Float kloc[])
    {
        this->fe = fe;
        this->shell = new ShellElement(fe->basis, fe->shape, thick, fe->xpts);
        int size = GetSize();
        for (auto npt = 0; npt < shell->n_qpts; npt++)
        {
            Float _kloc[size] = {0.0};
            shell->ComputeTransform(shell->qpts[npt]);
            this->operator()(_kloc);
            mat_ops::mat_mult_scalar(size, shell->J * shell->qWts[npt], _kloc);
            mat_ops::mat_add(size, true, _kloc, kloc, nullptr);
        }
        delete this->shell;
        // mat_ops::mat_print(15, 15, kloc);
        // exit(-1);
    }
}