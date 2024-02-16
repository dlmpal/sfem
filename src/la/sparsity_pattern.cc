#include "sparsity_pattern.h"
#include "petsc.h"

namespace sfem::la
{
    //=============================================================================
    std::pair<std::vector<int>, std::vector<int>> SparsityPattern(const mesh::Mesh *mesh, int n_vars)
    {
        // Pattern (node-to-node connectivity)
        struct _IntPairHash
        {
            size_t operator()(const std::pair<int, int> &p) const
            {
                auto hash1 = std::hash<int>{}(p.first);
                auto hash2 = std::hash<int>{}(p.second);
                return hash1 << 1 + hash1 ^ hash2;
            }
        };
        std::unordered_map<std::pair<int, int>, bool, _IntPairHash> pattern;

        // Number of non-zeros (nnz)
        std::vector<int> diag_nnz(mesh->GetNumNodesOwned(), 0);
        std::vector<int> off_diag_nnz(mesh->GetNumNodesOwned(), 0);
        std::vector<double> ghost_nnz(mesh->GetNumNodesGhost(), 0);
        for (const auto &cell : mesh->GetCells())
        {
            auto [nodes, _] = mesh->GetCellNodes(cell);
            for (auto i = 0; i < cell.n_nodes; i++)
            {
                for (auto j = 0; j < cell.n_nodes; j++)
                {
                    auto pair = std::make_pair(nodes[i], nodes[j]);
                    if (pattern.count(pair) == 0)
                    {
                        pattern[pair] = true;
                        // Node "i" is locally owned
                        if (mesh->IsNodeOwned(nodes[i]))
                        {
                            if (mesh->IsNodeOwned(nodes[j]))
                            {
                                diag_nnz[nodes[i]]++;
                            }
                            else
                            {
                                off_diag_nnz[nodes[i]]++;
                            }
                        }
                        // Node "i" is owned by another process (ghost)
                        else
                        {
                            ghost_nnz[nodes[i] - mesh->GetNumNodesOwned()]++;
                        }
                    }
                }
            }
        }

        // Add the diagonal non-zero contributions from other processes
        Vector ghost_nnz_vec(mesh->GetNumNodesOwned(), mesh->GetNumNodesGlobal(), mesh->GetGhostNodes(mesh::IndexingType::RENUMBERED));
        ghost_nnz_vec.AddValues(mesh->GetNumNodesGhost(), mesh->GetGhostNodes(mesh::IndexingType::RENUMBERED).data(), ghost_nnz.data());
        ghost_nnz_vec.Assemble();
        ghost_nnz = ghost_nnz_vec.GetLocalValues();
        for (auto i = 0; i < mesh->GetNumNodesOwned(); i++)
        {
            diag_nnz[i] += static_cast<int>(ghost_nnz[i]);
        }

        // Resize the NNZ vectors for the given number of variables per node
        diag_nnz.resize(mesh->GetNumNodesOwned() * n_vars);
        off_diag_nnz.resize(mesh->GetNumNodesOwned() * n_vars);
        for (auto i = mesh->GetNumNodesOwned() - 1; i >= 0; i--)
        {
            for (auto j = 0; j < n_vars; j++)
            {
                diag_nnz[i * n_vars + j] = diag_nnz[i] * n_vars;
                off_diag_nnz[i * n_vars + j] = off_diag_nnz[i] * n_vars;
            }
        }

        return std::make_pair(diag_nnz, off_diag_nnz);
    }
}