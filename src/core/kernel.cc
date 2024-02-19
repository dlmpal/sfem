#include "kernel.h"
#include "../common/mat_ops.h"

namespace sfem::kernel
{
    //=============================================================================
    Kernel::~Kernel(){};
    //=============================================================================
    int Kernel::GetSize() const
    {
        int n_dof = fe->GetNumDof();
        int size = n_dof * n_dof;
        if (type == KernelType::RHS)
        {
            size = n_dof;
        }
        else if (type == KernelType::BOTH)
        {
            size = n_dof * (n_dof + 1);
        }
        return size;
    }
    //=============================================================================
    std::vector<Float> Kernel::Integrate(fe::FiniteElement *fe)
    {
        this->fe = fe;
        int size = GetSize();
        auto basis = fe->GetBasis();
        std::vector<Float> kloc(size, 0.0);
        std::vector<Float> _kloc(size, 0.0);

        for (auto npt = 0; npt < basis->GetNumQuadPts(); npt++)
        {
            fe->ComputeTransform(basis->GetQuadPt(npt));
            Evaluate(_kloc);
            mat_ops::mat_mult_scalar(size, fe->J() * basis->GetQuadWeight(npt), _kloc.data());
            mat_ops::mat_add(size, true, _kloc.data(), kloc.data(), nullptr);
        }

        return kloc;
    }
}