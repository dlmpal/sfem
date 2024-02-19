#include "mass.h"

namespace sfem::kernel::basic
{
    //=============================================================================
    Mass::Mass(Float c)
        : c(c)
    {
    }
    //=============================================================================
    Kernel::KernelType Mass::GetType() const
    {
        return KernelType::MASS;
    }
    //=============================================================================
    Float Mass::GetCoefficient() const
    {
        return c;
    }
    //=============================================================================
    void Mass::Evaluate(std::vector<Float> &kloc)
    {
        auto n_vars = fe->GetNumVars();
        auto n_nodes = fe->GetNumNodes();
        auto n_dof = fe->GetNumDof();
        auto N = fe->N();

        for (auto i = 0; i < n_nodes; i++)
        {
            for (auto j = 0; j < n_vars; j++)
            {
                for (auto k = 0; k < n_nodes; k++)
                {
                    for (auto l = 0; l < n_vars; l++)
                    {
                        if (j == l)
                        {
                            kloc[(i * n_vars + j) * n_dof + (k * n_vars + l)] = c * N[i] * N[k];
                        }
                    }
                }
            }
        }
    }
}