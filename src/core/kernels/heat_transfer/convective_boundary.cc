#include "convective_boundary.h"

namespace sfem::kernel::heat_transfer
{
    ConvectiveBoundary::ConvectiveBoundary(Float h, Float T_inf, basic::Diffusion *kernel) : Kernel(ADD_BOTH)
    {
        this->h = h;
        this->T_inf = T_inf;
        this->kernel = kernel;
    }

    void ConvectiveBoundary::operator()(Float kloc[])
    {
        Float kappa = kernel->GetCoefficient();
        for (auto i = 0; i < fe->n_nodes; i++)
            for (auto j = 0; j < fe->n_nodes; j++)
                kloc[i * fe->n_nodes + j] += -kappa * h * fe->N[i] * fe->N[j];

        for (auto i = 0; i < fe->n_nodes; i++)
            kloc[fe->n_nodes * fe->n_nodes + i] = -kappa * h * T_inf * fe->N[i];
    }
}