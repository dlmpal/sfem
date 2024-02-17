#include "field.h"
#include "../common/error.h"

namespace sfem::field
{
    //=============================================================================
    Field::Field(const std::string &name, int n_vars, mesh::Mesh *mesh)
        : name(name), n_vars(n_vars), mesh(mesh)
    {
        if (mesh == nullptr)
        {
            /// @todo ERROR
        }

        n_dof = n_vars * mesh->GetNumNodes();
        n_dof_owned = n_vars * mesh->GetNumNodesOwned();
        n_dof_ghost = n_vars * mesh->GetNumNodesGhost();
        n_dof_global = n_vars * mesh->GetNumNodesGlobal();
        values.resize(n_dof);
    }
    //=============================================================================
    std::string Field::GetName() const
    {
        return name;
    }
    //=============================================================================
    int Field::GetNumVars() const
    {
        return n_vars;
    }
    //=============================================================================
    mesh::Mesh *Field::GetMesh() const
    {
        return mesh;
    }
    //=============================================================================
    int Field::GetNumDof() const
    {
        return n_dof;
    }
    //=============================================================================
    int Field::GetNumOwnedDof() const
    {
        return n_dof_owned;
    }
    //=============================================================================
    int Field::GetNumGhostDof() const
    {
        return n_dof_ghost;
    }
    //=============================================================================
    int Field::GetNumGlobalDof() const
    {
        return n_dof_global;
    }
    //=============================================================================
    void Field::GetNodeDof(int n_nodes, const int nodes[], int dof[]) const
    {
        for (auto i = 0; i < n_nodes; i++)
        {
            for (auto j = 0; j < n_vars; j++)
            {
                dof[i * n_vars + j] = nodes[i] * n_vars + j;
            }
        }
    }
    //=============================================================================
    std::vector<int> Field::GetOwnedDof(mesh::IndexingType indexing) const
    {
        auto owned_nodes = mesh->GetOwnedNodes(indexing);
        std::vector<int> owned_dof(n_dof_owned);
        GetNodeDof(owned_nodes.size(), owned_nodes.data(), owned_dof.data());
        return owned_dof;
    }
    //=============================================================================
    std::vector<int> Field::GetGhostDof(mesh::IndexingType indexing) const
    {
        auto ghost_nodes = mesh->GetGhostNodes(indexing);
        std::vector<int> ghost_dof(n_dof_ghost);
        GetNodeDof(ghost_nodes.size(), ghost_nodes.data(), ghost_dof.data());
        return ghost_dof;
    }
    //=============================================================================
    std::array<int, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD>
    Field::GetCellDof(const mesh::Cell &cell, mesh::IndexingType indexing) const
    {
        std::array<int, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD> cell_dof;
        auto [cell_nodes, _] = mesh->GetCellNodes(cell, indexing);
        GetNodeDof(cell.n_nodes, cell_nodes.data(), cell_dof.data());
        return cell_dof;
    }
    //=============================================================================
    std::array<Float, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD>
    Field::GetCellValues(const mesh::Cell &cell) const
    {
        std::array<Float, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD> cell_values;
        auto cell_dof = GetCellDof(cell, mesh::IndexingType::LOCAL);
        for (auto i = 0; i < cell.n_nodes * n_vars; i++)
        {
            cell_values[i] = values[cell_dof[i]];
        }

        return cell_values;
    }
    //=============================================================================
    void Field::AddFixedDof(const std::string &region_name, int var, Float value)
    {
        for (auto node : mesh->GetRegionNodes(region_name, mesh::IndexingType::RENUMBERED))
        {
            fixed_dof[node * n_vars + var] = value;
        }
    }
    //=============================================================================
    std::pair<std::vector<int>, std::vector<Float>> Field::GetFixedDof() const
    {
        std::vector<int> f_dof(fixed_dof.size());
        std::vector<Float> f_dof_values(fixed_dof.size());
        for (auto kv : fixed_dof)
        {
            f_dof.push_back(kv.first);
            f_dof_values.push_back(kv.second);
        }
        return std::make_pair(f_dof, f_dof_values);
    }
    //=============================================================================
    void Field::ClearFixedDof()
    {
        fixed_dof.clear();
    }
    //=============================================================================
    void Field::SetValues(const std::vector<Float> &values)
    {
        if (values.size() != n_dof)
        {
            error::InvalidSizeError(n_dof, values.size(), __FILE__, __LINE__);
        }
        for (auto i = 0; i < n_dof; i++)
        {
            this->values[i] = values[i];
        }
    }
    //=============================================================================
    const std::vector<Float> &Field::GetValues() const
    {
        return values;
    }
    //=============================================================================
    std::vector<Float> AssembleGlobalValues(const Field &field)
    {
        int n_procs = Logger::GetInstance().GetNumProcs();
        int proc_rank = Logger::GetInstance().GetProcRank();

        // If serial, return early
        if (n_procs == 1)
        {
            return field.GetValues();
        }
#ifdef SFEM_USE_MPI
        else
        {
            // Get number of dof per process
            std::vector<int> n_dof_per_proc(n_procs);
            int n_dof_owned = field.GetNumOwnedDof();
            MPI_Gather(&n_dof_owned, 1, MPI_INT, n_dof_per_proc.data(),
                       1, MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);

            // Get the dof offsets
            std::vector<int> dof_offsets(n_procs);
            int dof_offset = 0;
            for (auto i = 0; i < n_procs; i++)
            {
                dof_offsets[i] = dof_offset;
                dof_offset += n_dof_per_proc[i];
            }

            // Get the dof owned by this process
            auto owned_dof = field.GetOwnedDof(mesh::IndexingType::GLOBAL);

            // Assemble the owned dof in one big vector
            std::vector<int> all_dof;
            if (Logger::GetInstance().GetProcRank() == SFEM_ROOT)
            {
                all_dof.resize(field.GetNumGlobalDof());
            }
            MPI_Gatherv(owned_dof.data(), field.GetNumOwnedDof(), MPI_INT, all_dof.data(), n_dof_per_proc.data(),
                        dof_offsets.data(), MPI_INT, SFEM_ROOT, SFEM_COMM_WORLD);

            // Finally, assemble all the values in one vector, and re-order it
            std::vector<Float> all_values;
            if (proc_rank == SFEM_ROOT)
            {
                all_values.resize(field.GetNumGlobalDof());
            }
            MPI_Gatherv(field.GetValues().data(), field.GetNumOwnedDof(), SFEM_MPI_FLOAT, all_values.data(),
                        n_dof_per_proc.data(), dof_offsets.data(), SFEM_MPI_FLOAT, SFEM_ROOT, SFEM_COMM_WORLD);

            std::vector<Float> all_values_reordered;
            if (proc_rank == SFEM_ROOT)
            {
                all_values_reordered.resize(field.GetNumGlobalDof());
                for (auto i = 0; i < field.GetNumGlobalDof(); i++)
                {
                    all_values_reordered[all_dof[i]] = all_values[i];
                }
            }
            return all_values_reordered;
        }
#endif // SFEM_USE_MPI
    }

    // void FieldGrad(field::Field *field)
    // {
    //     std::vector<field::Field *> grad(field->n_vars);
    //     for (auto i = 0; i < field->n_vars; i++)
    //         grad[i] = field::CreateField("grad_" + field->name + "_" + std::to_string(i),
    //                                      field->mesh->dim,
    //                                      field->mesh);

    //     la::PetscVec contribs_per_node((field->n_dof_owned + field->n_dof_ghost) / (field->n_vars),
    //                                    (field->n_dof) / (field->n_vars),
    //                                    (field->n_dof_ghost) / (field->n_vars),
    //                                    field->GetGhostDof().data());

    //     for (auto cell : field->mesh->cells)
    //     {
    //         auto [nodes, xpts] = field->mesh->GetCellNodes(cell, mesh::Mesh::RENUMBERED_INDEXING);

    //         std::vector<std::array<int, 300>> dof(grad.size());
    //         for (auto i = 0; i < field->n_vars; i++)
    //             dof[i] = grad[i]->GetCellDof(cell);

    //         auto values = field->GetCellValues(cell);
    //         fe::FiniteElement fe(cell, field->mesh->dim, field->n_vars, xpts.data());

    //         if (fe.basis->dim < field->mesh->dim)
    //             continue;

    //         for (auto i = 0; i < fe.n_nodes; i++)
    //         {
    //             Float cell_grad[fe.n_vars * fe.dim] = {0};
    //             fe.ComputeFieldGrad(fe.basis->xi[i], values.data(), cell_grad);
    //             for (auto j = 0; j < fe.n_vars; j++)
    //                 for (auto k = 0; k < fe.dim; k++)
    //                     grad[j]->values[dof[j][i * grad[j]->n_vars + k]] += cell_grad[j * fe.dim + k];
    //         }

    //         Float ones[cell.n_nodes];
    //         mat_ops::mat_set(cell.n_nodes, 1.0, ones);
    //         contribs_per_node.AddValues(cell.n_nodes, nodes.data(), ones);
    //     }

    //     for (auto j = 0; j < field->n_vars; j++)
    //         for (auto i = 0; i < field->mesh->n_nodes; i++)
    //             for (auto k = 0; k < field->mesh->dim; k++)
    //                 grad[j]->values[i * field->mesh->dim + k] = grad[j]->values[i * field->mesh->dim + k] / contribs_per_node[i];

    //     for (auto i = 0; i < field->n_vars; i++)
    //         grad[i]->WriteValues();
    // }
}