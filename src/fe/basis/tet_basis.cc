#include "basis.h"

namespace sfem::fe::basis
{
    //=============================================================================
    Tet1Basis::Tet1Basis()
    {
        dim = 3;
        n_nodes = 4;
        n_qpts = 1;

        // Corner vertices
        xi[0][0] = 0.0;
        xi[0][1] = 0.0;
        xi[0][2] = 0.0;

        xi[1][0] = 1.0;
        xi[1][1] = 0.0;
        xi[1][2] = 0.0;

        xi[2][0] = 0.0;
        xi[2][1] = 1.0;
        xi[2][2] = 0.0;

        xi[3][0] = 0.0;
        xi[3][1] = 0.0;
        xi[3][2] = 0.1;

        // Quadrature point
        qpts[0][0] = 0.25;
        qpts[0][1] = 0.25;
        qpts[0][2] = 0.25;

        // Quadrature weight
        qwts[0] = 0.166666667;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    Tet1Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        N[0] = 1 - pt[0] - pt[1] - pt[2];
        N[1] = pt[0];
        N[2] = pt[1];
        N[3] = pt[2];

        dNdxi[0][0] = -1.0;
        dNdxi[0][1] = -1.0;
        dNdxi[0][2] = -1.0;

        dNdxi[1][0] = 1.0;
        dNdxi[1][1] = 0.0;
        dNdxi[1][2] = 0.0;

        dNdxi[2][0] = 0.0;
        dNdxi[2][1] = 1.0;
        dNdxi[2][2] = 0.0;

        dNdxi[3][0] = 0.0;
        dNdxi[3][1] = 0.0;
        dNdxi[3][2] = 1.0;

        return std::make_pair(N, dNdxi);
    };
    //=============================================================================
    Tet2Basis::Tet2Basis()
    {
        dim = 3;
        n_nodes = 10;
        n_qpts = 5;

        // Corner vertices
        xi[0][0] = 0.0;
        xi[0][1] = 0.0;
        xi[0][2] = 0.0;

        xi[1][0] = 1.0;
        xi[1][1] = 0.0;
        xi[1][2] = 0.0;

        xi[2][0] = 0.0;
        xi[2][1] = 1.0;
        xi[2][2] = 0.0;

        xi[3][0] = 0.0;
        xi[3][1] = 0.0;
        xi[3][2] = 0.1;

        // Mid-side vertices
        xi[4][0] = 0.5;
        xi[4][1] = 0.0;
        xi[4][2] = 0.0;

        xi[5][0] = 0.5;
        xi[5][1] = 0.5;
        xi[5][2] = 0.0;

        xi[6][0] = 0.0;
        xi[6][1] = 0.5;
        xi[6][2] = 0.0;

        xi[7][0] = 0.0;
        xi[7][1] = 0.0;
        xi[7][2] = 0.5;

        xi[8][0] = 0.0;
        xi[8][1] = 0.5;
        xi[8][2] = 0.5;

        xi[9][0] = 0.5;
        xi[9][1] = 0.0;
        xi[9][2] = 0.5;

        // Quadrature points
        qpts[0][0] = 0.25;
        qpts[0][1] = 0.25;
        qpts[0][2] = 0.25;

        qpts[1][0] = 1.0 / 6.0;
        qpts[1][1] = 1.0 / 6.0;
        qpts[1][2] = 1.0 / 6.0;

        qpts[2][0] = 0.5;
        qpts[2][1] = 1.0 / 6.0;
        qpts[2][2] = 1.0 / 6.0;

        qpts[3][0] = 1.0 / 6.0;
        qpts[3][1] = 0.5;
        qpts[3][2] = 1.0 / 6.0;

        qpts[4][0] = 1.0 / 6.0;
        qpts[4][1] = 1.0 / 6.0;
        qpts[4][2] = 0.5;

        // Quadrature Weights
        qwts[0] = -2.0 / 15.0;
        qwts[1] = 3.0 / 40.0;
        qwts[2] = 3.0 / 40.0;
        qwts[3] = 3.0 / 40.0;
        qwts[4] = 3.0 / 40.0;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    Tet2Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        Float l0 = 1 - pt[0] - pt[1] - pt[2];
        Float l1 = pt[0];
        Float l2 = pt[1];
        Float l3 = pt[2];

        // Corner vertices
        N[0] = l0 * (2 * l0 - 1);
        N[1] = l1 * (2 * l1 - 1);
        N[2] = l2 * (2 * l2 - 1);
        N[3] = l3 * (2 * l3 - 1);

        // Mid-side vertices
        N[4] = 4 * l1 * l0;
        N[5] = 4 * l1 * l2;
        N[6] = 4 * l2 * l0;
        N[7] = 4 * l3 * l0;
        N[8] = 4 * l2 * l3;
        N[9] = 4 * l1 * l3;

        // Corner vertices
        dNdxi[0][0] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
        dNdxi[0][1] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;
        dNdxi[0][2] = 4 * pt[0] + 4 * pt[1] + 4 * pt[2] - 3;

        dNdxi[1][0] = 4 * pt[0] - 1;
        dNdxi[1][1] = 0;
        dNdxi[1][2] = 0;

        dNdxi[2][0] = 0;
        dNdxi[2][1] = 4 * pt[1] - 1;
        dNdxi[2][2] = 0;

        dNdxi[3][0] = 0;
        dNdxi[3][1] = 0;
        dNdxi[3][2] = 4 * pt[2] - 1;

        // Mid-side vertices
        dNdxi[4][0] = -4 * (2 * pt[0] + pt[1] + pt[2] - 1);
        dNdxi[4][1] = -4 * pt[0];
        dNdxi[4][2] = -4 * pt[0];

        dNdxi[5][0] = 4 * pt[1];
        dNdxi[5][1] = 4 * pt[0];
        dNdxi[5][2] = 0.0;

        dNdxi[6][0] = -4 * pt[1];
        dNdxi[6][1] = -4 * (pt[0] + 2 * pt[1] + pt[2] - 1);
        dNdxi[6][2] = -4 * pt[1];

        dNdxi[7][0] = -4 * pt[2];
        dNdxi[7][1] = -4 * pt[2];
        dNdxi[7][2] = -4 * (pt[0] + pt[1] + 2 * pt[2] - 1);

        dNdxi[8][0] = 0;
        dNdxi[8][1] = 4 * pt[2];
        dNdxi[8][2] = 4 * pt[1];

        dNdxi[9][0] = 4 * pt[2];
        dNdxi[9][1] = 0;
        dNdxi[9][2] = 4 * pt[0];

        return std::make_pair(N, dNdxi);
    };

}