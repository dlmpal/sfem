#include "pressure.h"

namespace sfem::kernel::elasticity
{
    //=============================================================================
    Pressure::Pressure(field::Field *p) : Kernel(KernelType::RHS)
    {
        this->p = p;
    };
    //=============================================================================
    void Pressure::operator()(Float kloc[])
    {
        auto values = p->GetCellValues(fe->GetCell());
        auto normal = fe->GetShape()->FaceNormal(-1, fe->GetXpts().data());
        auto n_vars = fe->GetNumVars();
        auto n_nodes = fe->GetNumNodes();
        auto N = fe->N();
        for (auto i = 0; i < n_nodes; i++)
            for (auto j = 0; j < n_vars; j++)
                kloc[i * n_vars + j] = -values[i] * normal.x[j] * N[i];
    }
}