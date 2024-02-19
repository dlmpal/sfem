#include "diffusion.h"

namespace sfem::kernel::basic
{
    //=============================================================================
    Diffusion::Diffusion(Float c)
        : c(c)
    {
    }
    //=============================================================================
    Kernel::KernelType Diffusion::GetType() const
    {
        return KernelType::LHS;
    }
    //=============================================================================
    Float Diffusion::GetCoefficient() const
    {
        return c;
    }
    //=============================================================================
    void Diffusion::Evaluate(std::vector<Float> &kloc)
    {
        auto n_nodes = fe->GetNumNodes();
        auto dNdX = fe->dNdX();

        for (int i = 0; i < n_nodes; i++)
        {
            for (int j = 0; j < n_nodes; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    kloc[i * n_nodes + j] += c * dNdX[i][k] * dNdX[j][k];
                }
            }
        }
    }
}