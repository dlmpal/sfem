#include "basis.h"

namespace sfem::fe::basis
{
    //=============================================================================
    L1Basis::L1Basis()
    {
        dim = 1;
        n_nodes = 2;
        n_qpts = 1;

        xi[0][0] = -1;
        xi[1][0] = 1;

        qpts[0][0] = 0.0;

        qwts[0] = 2.0;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    L1Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        N[0] = 0.5 * (1.0 - pt[0]);
        N[1] = 0.5 * (1.0 + pt[0]);

        dNdxi[0][0] = -0.5;
        dNdxi[1][0] = 0.5;

        return std::make_pair(N, dNdxi);
    }
    //=============================================================================
    L2Basis::L2Basis()
    {
        dim = 1;
        n_nodes = 3;
        n_qpts = 2;

        xi[0][0] = -1;
        xi[1][0] = 0;
        xi[2][0] = 1;

        qpts[0][0] = -0.577350269189626;
        qpts[0][1] = 0.577350269189626;

        qwts[0] = 1.0;
        qwts[1] = 1.0;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    L2Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        N[0] = -0.5 * pt[0] * (1.0 - pt[0]);
        N[1] = (1.0 - pt[0]) * (1.0 + pt[0]);
        N[2] = 0.5 * (1.0 + pt[0]) * pt[0];

        dNdxi[0][0] = -0.5 + pt[0];
        dNdxi[1][0] = -2.0 * pt[0];
        dNdxi[2][0] = 0.5 + pt[0];

        return std::make_pair(N, dNdxi);
    }

}