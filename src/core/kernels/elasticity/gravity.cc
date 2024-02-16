#include "gravity.h"

namespace sfem::kernel::elasticity
{
    //=============================================================================
    Gravity::Gravity(Float rho, int direction) : Kernel(KernelType::RHS)
    {
        this->rho = rho;
        this->direction = direction;
    }
    //=============================================================================
    void Gravity::operator()(Float kloc[])
    {
        auto n_vars = fe->GetNumVars();
        auto n_nodes = fe->GetNumNodes();
        auto N = fe->N();
        for (auto i = 0; i < n_nodes; i++)
            kloc[i * n_vars + direction] = -rho * g * N[i];
    }
}