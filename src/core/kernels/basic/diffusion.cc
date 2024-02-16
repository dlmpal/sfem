#include "diffusion.h"

namespace sfem::kernel::basic
{
    //=============================================================================
    Diffusion::Diffusion(Float D) : Kernel(KernelType::LHS)
    {
        c = D;
    }
    //=============================================================================
    Float Diffusion::GetCoefficient() const
    {
        return c;
    }
    //=============================================================================
    void Diffusion::operator()(Float kloc[])
    {
        auto n_nodes = fe->GetNumNodes();
        auto dNdX = fe->dNdX();
        for (int i = 0; i < n_nodes; i++)
            for (int j = 0; j < n_nodes; j++)
                for (int k = 0; k < 3; k++)
                    kloc[i * n_nodes + j] += c * dNdX[i][k] * dNdX[j][k];
    }
}