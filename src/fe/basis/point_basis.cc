#include "basis.h"

namespace sfem::fe::basis
{
    //=============================================================================
    PointBasis::PointBasis()
    {
        dim = 0;
        n_nodes = 1;
        n_qpts = 1;

        xi[0][0] = 1.0;

        qpts[0][0] = 1.0;

        qwts[0] = 1.0;
    }
    //=============================================================================
    std::pair<std::array<Float, SFEM_MAX_CELL_NODES>, std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES>>
    PointBasis::ComputeShapeAndGrad(const std::array<Float, 3> &pt) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES> N;
        std::array<std::array<Float, 3>, SFEM_MAX_CELL_NODES> dNdxi;

        N[0] = 1.0;
        dNdxi[0][0] = 1.0;

        return std::make_pair(N, dNdxi);
    }

}
