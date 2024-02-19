#include "source.h"

namespace sfem::kernel::basic
{
    //=============================================================================
    Source::Source(const Float c[])
        : c(c)
    {
    }
    //=============================================================================
    Kernel::KernelType Source::GetType() const
    {
        return KernelType::RHS;
    }
    //=============================================================================
    void Source::ComputeValues(Float values[])
    {
        auto n_vars = fe->GetNumVars();
        auto n_nodes = fe->GetNumNodes();
        for (auto i = 0; i < n_nodes; i++)
        {
            for (auto j = 0; j < n_vars; j++)
            {
                values[i * n_vars + j] = c[j];
            }
        }
    }
    //=============================================================================
    void Source::Evaluate(std::vector<Float> &kloc)
    {
        auto n_vars = fe->GetNumVars();
        auto n_nodes = fe->GetNumNodes();
        auto n_dof = fe->GetNumDof();
        Float values[n_dof];
        ComputeValues(values);
        auto N = fe->N();

        for (auto i = 0; i < n_nodes; i++)
        {
            for (auto j = 0; j < n_vars; j++)
            {
                kloc[i * n_vars + j] = values[i * n_vars + j] * N[i];
            }
        }
    }
}