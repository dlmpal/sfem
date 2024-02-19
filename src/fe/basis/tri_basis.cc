#include "basis.h"

namespace sfem::fe::basis
{
    //=============================================================================
    Tri1Basis::Tri1Basis()
    {
        dim = 2;
        n_nodes = 3;
        n_qpts = 1;

        // Corner vertices
        xi[0][0] = 0.0;
        xi[0][1] = 0.0;

        xi[1][0] = 1.0;
        xi[1][1] = 0.0;

        xi[2][0] = 0.0;
        xi[2][1] = 1.0;

        // Quadrature point
        qpts[0][0] = 1.0 / 3.0;
        qpts[0][1] = 1.0 / 3.0;

        // Quadrature weight
        qwts[0] = 1.0;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    Tri1Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;
        N[0] = 1 - pt[0] - pt[1];
        N[1] = pt[0];
        N[2] = pt[1];

        dNdxi[0][0] = -1.0;
        dNdxi[0][1] = -1.0;

        dNdxi[1][0] = 1.0;
        dNdxi[1][1] = 0.0;

        dNdxi[2][0] = 0.0;
        dNdxi[2][1] = 1.0;

        return std::make_pair(N, dNdxi);
    };
    //=============================================================================
    Tri2Basis::Tri2Basis()
    {
        dim = 2;
        n_nodes = 6;
        n_qpts = 6;

        // Corner vertices
        xi[0][0] = 0.0;
        xi[0][1] = 0.0;

        xi[1][0] = 1.0;
        xi[1][1] = 0.0;

        xi[2][0] = 0.0;
        xi[2][1] = 1.0;

        // Mid-point vertices
        xi[3][0] = 0.5;
        xi[3][1] = 0.0;

        xi[4][0] = 0.5;
        xi[4][1] = 0.5;

        xi[5][0] = 0.0;
        xi[5][1] = 0.5;

        // Quadrature points
        qpts[0][0] = 0.091576213509771;
        qpts[0][1] = 0.091576213509771;

        qpts[1][0] = 0.816847572980459;
        qpts[1][1] = 0.091576213509771;

        qpts[2][0] = 0.091576213509771;
        qpts[2][1] = 0.816847572980459;

        qpts[3][0] = 0.108103018168070;
        qpts[3][1] = 0.108103018168070;

        qpts[4][0] = 0.445948490915965;
        qpts[4][1] = 0.108103018168070;

        qpts[5][0] = 0.108103018168070;
        qpts[5][1] = 0.445948490915965;

        // Quadrature weights
        qwts[0] = 0.054975871827661;
        qwts[1] = 0.054975871827661;
        qwts[2] = 0.054975871827661;
        qwts[3] = 0.111690794839006;
        qwts[4] = 0.111690794839006;
        qwts[5] = 0.111690794839006;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    Tri2Basis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        const Float l0 = 1 - pt[0] - pt[1];
        const Float l1 = pt[0];
        const Float l2 = pt[1];

        N[0] = l0 * (2 * l0 - 1);
        N[1] = l1 * (2 * l1 - 1);
        N[2] = l2 * (2 * l2 - 1);
        N[3] = 4 * l0 * l1;
        N[4] = 4 * l1 * l2;
        N[5] = 4 * l0 * l2;

        dNdxi[0][0] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;
        dNdxi[0][1] = 4.0 * pt[0] + 4.0 * pt[1] - 3.0;

        dNdxi[1][0] = 4.0 * pt[0] - 1.0;
        dNdxi[1][1] = 0.0;

        dNdxi[2][0] = 0.0;
        dNdxi[2][1] = 4.0 * pt[1] - 1.0;

        dNdxi[3][0] = 4.0 - 8.0 * pt[0] - 4.0 * pt[1];
        dNdxi[3][1] = -4.0 * pt[0];

        dNdxi[4][0] = 4.0 * pt[1];
        dNdxi[4][1] = 4.0 * pt[0];

        dNdxi[5][0] = -4.0 * pt[1];
        dNdxi[5][1] = 4.0 - 4.0 * pt[0] - 8.0 * pt[1];

        return std::make_pair(N, dNdxi);
    };
}