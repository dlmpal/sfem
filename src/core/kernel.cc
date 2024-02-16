#include "kernel.h"
#include "../common/mat_ops.h"

namespace sfem::kernel
{
    //=============================================================================
    Kernel::Kernel(KernelType type)
        : type(type)
    {
    }
    //=============================================================================
    Kernel::~Kernel(){

    };
    //=============================================================================
    Kernel::KernelType Kernel::GetType() const
    {
        return type;
    }
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
    void Kernel::Integrate(fe::FiniteElement *fe, Float kloc[])
    {
        this->fe = fe;
        int size = GetSize();
        auto basis = fe->GetBasis();
        for (auto npt = 0; npt < basis->GetNumQuadPts(); npt++)
        {
            Float _kloc[size] = {0.0};
            fe->ComputeTransform(basis->GetQuadPt(npt));
            this->operator()(_kloc);
            mat_ops::mat_mult_scalar(size, fe->J() * basis->GetQuadWeight(npt), _kloc);
            mat_ops::mat_add(size, true, _kloc, kloc, nullptr);
        }
    }

}