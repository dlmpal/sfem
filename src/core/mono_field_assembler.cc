#include "mono_field_assembler.h"
#include "../common/mat_ops.h"
#include "chrono"
#include "../common/logger.h"

namespace sfem::assembly
{
    //=============================================================================
    MonoFieldAssembler::MonoFieldAssembler(field::Field *field, AssemblyType type)
        : field(field), type(type)
    {
        if (field == nullptr)
        {
            // ERROR
        }

        auto [diag_nnz, offdiag_nnz] = la::SparsityPattern(field->GetMesh(), field->GetNumVars());

        U = std::unique_ptr<la::Vector>(new la::Vector(field->GetNumOwnedDof(), field->GetNumGlobalDof(), field->GetGhostDof()));
        F = std::unique_ptr<la::Vector>(new la::Vector(field->GetNumOwnedDof(), field->GetNumGlobalDof(), field->GetGhostDof()));
        K = std::unique_ptr<la::Matrix>(new la::Matrix(field->GetNumOwnedDof(), field->GetNumGlobalDof(), diag_nnz, offdiag_nnz));
        M = nullptr;
        if (type == AssemblyType::DYNAMIC)
        {
            M = std::unique_ptr<la::Matrix>(new la::Matrix(field->GetNumOwnedDof(), field->GetNumGlobalDof(), diag_nnz, offdiag_nnz));
        }
    }
    //=============================================================================
    AssemblyType MonoFieldAssembler::GetType() const
    {
        return type;
    }
    //=============================================================================
    field::Field *MonoFieldAssembler::GetField() const
    {
        return field;
    }
    //=============================================================================
    la::Vector *MonoFieldAssembler::GetSolutionVector() const
    {
        return U.get();
    }
    //=============================================================================
    la::Vector *MonoFieldAssembler::GetLoadVector() const
    {
        return F.get();
    }
    //=============================================================================
    la::Matrix *MonoFieldAssembler::GetStiffnessMatrix() const
    {
        return K.get();
    }
    //=============================================================================
    la::Matrix *MonoFieldAssembler::GetMassMatrix() const
    {
        return M.get();
    }
    //=============================================================================
    void MonoFieldAssembler::AddKernel(std::string region, kernel::Kernel *kernel)
    {
        kernel_map[region].push_back(kernel);
    }
    //=============================================================================
    void MonoFieldAssembler::AssembleSystem()
    {
        // Reset vectors and matrices
        F->SetAll(0.0);
        U->SetAll(0.0);
        K->Reset();
        if (type == AssemblyType::DYNAMIC)
        {
            M->Reset();
        }

        // Update solution vector values
        const auto &values = field->GetValues();
        U->InsertValues(field->GetNumOwnedDof(), field->GetOwnedDof().data(), values.data());

        auto mesh = field->GetMesh();
        for (auto region : mesh->GetRegions())
        {
            // If no kernels are defined for this region, continue
            if (kernel_map.count(region.name) == 0)
            {
                continue;
            }

            // Integrate all kernels for this region
            for (const auto &cell : mesh->GetRegionCells(region.name))
            {
                auto [_, xpts] = mesh->GetCellNodes(cell);
                auto dof = field->GetCellDof(cell);
                fe::FiniteElement fe(cell, mesh->GetDim(), field->GetNumVars(), xpts);

                for (auto kernel : kernel_map[region.name])
                {
                    auto kloc = kernel->Integrate(&fe);

                    switch (kernel->GetType())
                    {
                    case kernel::Kernel::KernelType::LHS:
                        K->AddValues(fe.GetNumDof(), dof.data(), kloc.data());
                        break;
                    case kernel::Kernel::KernelType::MASS:
                        if (type == AssemblyType::DYNAMIC)
                            M->AddValues(fe.GetNumDof(), dof.data(), kloc.data());
                        break;
                    case kernel::Kernel::KernelType::RHS:
                        F->AddValues(fe.GetNumDof(), dof.data(), kloc.data());
                        break;
                    case kernel::Kernel::KernelType::BOTH:
                        K->AddValues(fe.GetNumDof(), dof.data(), kloc.data());
                        break;
                    }
                }
            }
        }

        // Eliminate the fixed DOF and assemble the linear system
        auto fixed_dof = field->GetFixedDof();
        K->Assemble();
        K->ZeroRowsColumns(fixed_dof.first.size(), fixed_dof.first.data(), fixed_dof.second.data(), U.get(), F.get());
        U->Assemble();
        F->Assemble();
        if (type == AssemblyType::DYNAMIC)
        {
            M->Assemble();
        }
    };
    //=============================================================================
    void MonoFieldAssembler::UpdateFieldValues()
    {
        field->SetValues(U->GetLocalValues());
    }
}