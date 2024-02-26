#include "mesh_partitioner.h"
#include "../common/logger.h"
#include "../common/error.h"
#include "../third_party/sfem_metis.h"

namespace sfem::mesh
{
    //=============================================================================
    MeshPartitioner::MeshPartitioner(const Mesh *mesh, int n_parts)
        : mesh(mesh), n_parts(n_parts)
    {
        if (mesh == nullptr)
        {
            // ERROR
        }
    }
    //=============================================================================
    void MeshPartitioner::PartitionMesh()
    {
        auto [cpart, npart] = CreatePartition();
        AssignCells(cpart);
        AssignNodes(npart);
        RenumberNodes();
    }
    //=============================================================================
    void MeshPartitioner::AssignCells(std::vector<int> cpart)
    {
        n_cell_node_idx_per_part = std::vector<int>(n_parts, 0);

        // Count cells per partition
        int n_cells_per_part[n_parts] = {0};
        const auto &cells = mesh->GetCells();
        for (auto i = 0; i < mesh->GetNumCells(); i++)
        {
            n_cells_per_part[cpart[i]]++;
            n_cell_node_idx_per_part[cpart[i]] += cells[i].n_nodes;
        }

        // Assign cells to each partition
        cells_per_part.clear();
        for (auto i = 0; i < n_parts; i++)
        {
            cells_per_part.push_back(std::vector<int>());
            cells_per_part[i].reserve(n_cells_per_part[i]);
        }
        for (auto i = 0; i < mesh->GetNumCells(); i++)
        {
            cells_per_part[cpart[i]].push_back(i);
        }
    }
    //=============================================================================
    void MeshPartitioner::AssignNodes(std::vector<int> npart)
    {
        // Count owned nodes per partition
        int n_nodes_owned_per_part[n_parts] = {0};
        for (auto i = 0; i < mesh->GetNumNodes(); i++)
        {
            n_nodes_owned_per_part[npart[i]]++;
        }

        // Assign owned nodes
        nodes_owned_per_part.clear();
        for (auto i = 0; i < n_parts; i++)
        {
            nodes_owned_per_part.push_back(std::vector<int>());
            nodes_owned_per_part[i].reserve(n_nodes_owned_per_part[i]);
        }
        for (auto i = 0; i < mesh->GetNumNodes(); i++)
        {
            nodes_owned_per_part[npart[i]].push_back(i);
        }

        // Count ghost nodes per partition
        const auto &cells = mesh->GetCells();
        const auto &cell_node_conn = mesh->GetCellNodeConn();
        int n_nodes_ghost_per_proc[n_parts] = {0};
        for (auto i = 0; i < n_parts; i++)
        {
            std::unordered_map<int, short> is_node_included;
            for (auto j = 0; j < cells_per_part[i].size(); j++)
            {
                Cell cell = cells[cells_per_part[i][j]];
                for (auto k = 0; k < cell.n_nodes; k++)
                {
                    int node = cell_node_conn[cell.first_node_idx + k];
                    if (is_node_included.count(node) == 0 && npart[node] != i)
                    {
                        is_node_included[node] = 1;
                        n_nodes_ghost_per_proc[i]++;
                    }
                }
            }
        }

        // Assign ghost nodes
        nodes_ghost_per_part.clear();
        for (auto i = 0; i < n_parts; i++)
        {
            nodes_ghost_per_part.push_back(std::vector<int>());
            nodes_ghost_per_part[i].reserve(n_nodes_ghost_per_proc[i]);
        }
        for (auto i = 0; i < n_parts; i++)
        {
            std::unordered_map<int, short> is_node_included;
            for (auto j = 0; j < cells_per_part[i].size(); j++)
            {
                auto cell = cells[cells_per_part[i][j]];
                for (auto k = 0; k < cell.n_nodes; k++)
                {
                    int node = cell_node_conn[cell.first_node_idx + k];
                    if (is_node_included.count(node) == 0 && npart[node] != i)
                    {
                        is_node_included[node] = 1;
                        nodes_ghost_per_part[i].push_back(node);
                    }
                }
            }
        }
    }
    //=============================================================================
    void MeshPartitioner::RenumberNodes()
    {
        // Renumber owned nodes
        int node_offsets[n_parts], offset = 0;
        nodes_owned_per_part_renumbered.clear();
        for (auto i = 0; i < n_parts; i++)
        {
            node_offsets[i] = offset;
            offset += nodes_owned_per_part[i].size();
            nodes_owned_per_part_renumbered.push_back(std::vector<int>(nodes_owned_per_part[i].size()));
        }
        std::vector<int> renumbered(mesh->GetNumNodes());
        for (auto i = 0; i < n_parts; i++)
        {
            for (auto j = 0; j < nodes_owned_per_part[i].size(); j++)
            {
                renumbered[nodes_owned_per_part[i][j]] = node_offsets[i] + j;
                nodes_owned_per_part_renumbered[i][j] = node_offsets[i] + j;
            }
        }

        // Renumber ghost nodes
        nodes_ghost_per_part_renumbered.clear();
        for (auto i = 0; i < n_parts; i++)
        {
            nodes_ghost_per_part_renumbered.push_back(std::vector<int>(nodes_ghost_per_part[i].size()));
        }
        for (auto i = 0; i < n_parts; i++)
        {
            for (auto j = 0; j < nodes_ghost_per_part[i].size(); j++)
            {
                nodes_ghost_per_part_renumbered[i][j] = renumbered[nodes_ghost_per_part[i][j]];
            }
        }
    }
    //=============================================================================
    void MeshPartitioner::WritePartition(const std::string &path)
    {
        WriteCellPartition(path);
        WriteNodePartition(path);
    }
    //=============================================================================
    void MeshPartitioner::WriteCellPartition(const std::string &path)
    {
        std::ofstream file(path + "/CellPartition");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/CellPartition", __FILE__, __LINE__);
        }

        // File header
        file << n_parts << "\n";
        for (auto i = 0; i < n_parts; i++)
        {
            file << cells_per_part[i].size() << " " << n_cell_node_idx_per_part[i] << "\n";
        }

        // File body
        for (auto i = 0; i < n_parts; i++)
        {
            for (auto j = 0; j < cells_per_part[i].size(); j++)
            {
                file << cells_per_part[i][j] << "\n";
            }
        }
    }
    //=============================================================================
    void MeshPartitioner::WriteNodePartition(const std::string &path)
    {
        std::ofstream file(path + "/NodePartition");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/NodePartition", __FILE__, __LINE__);
        }

        // File header
        file << n_parts << "\n";
        for (auto i = 0; i < n_parts; i++)
        {
            file << nodes_owned_per_part[i].size() << " " << nodes_ghost_per_part[i].size() << "\n";
        }

        // File body
        for (auto i = 0; i < n_parts; i++)
        {
            for (auto j = 0; j < nodes_owned_per_part[i].size(); j++)
            {
                file << nodes_owned_per_part[i][j] << " " << nodes_owned_per_part_renumbered[i][j] << "\n";
            }
            for (auto j = 0; j < nodes_ghost_per_part[i].size(); j++)
            {
                file << nodes_ghost_per_part[i][j] << " " << nodes_ghost_per_part_renumbered[i][j] << "\n";
            }
        }
    }
    //=============================================================================
    MeshPartitioner *CreateMeshPartitioner(const std::string &type, const mesh::Mesh *mesh, int n_parts)
    {
        MeshPartitioner *partitioner = nullptr;

        if (type == "METIS")
        {
#ifdef SFEM_USE_METIS
            partitioner = new METISPartitioner(mesh, n_parts);
#else
            Logger::GetInstance().Error("SFEM was not compiled with METIS. Add SFEM_USE_METIS=ON and re-compile the library.\n", __FILE__, __LINE__);
#endif // SFEM_USE_METIS
        }

        if (partitioner == nullptr)
        {
            // ERROR
        }

        return partitioner;
    }
    //=============================================================================
#ifdef SFEM_USE_METIS
    METISPartitioner::METISPartitioner(const Mesh *mesh, int n_parts) : MeshPartitioner(mesh, n_parts)
    {
    }
    //=============================================================================
    std::pair<std::vector<int>, std::vector<int>> METISPartitioner::CreatePartition()
    {
        // METIS-required arrays
        // The eptr array contains the first and last node of each cell
        // The eind array contains the nodes corresponding to each cell
        // For example for element i eind[eptr[i] + j] is the node index of the j-th node for cell i
        // The cpart and npart arrays contain the cell and node partitions respectively
        auto n_cells = mesh->GetNumCells(), n_nodes = mesh->GetNumNodes(), n_cell_node_idx = mesh->GetSizeCellNodeConn();
        std::vector<idx_t> eptr(n_cells + 1);
        std::vector<idx_t> eind(n_cell_node_idx);
        std::vector<idx_t> cpart(n_cells);
        std::vector<idx_t> npart(n_nodes);

        const auto &cells = mesh->GetCells();
        const auto &cell_node_conn = mesh->GetCellNodeConn();
        for (auto i = 0; i < n_cells; i++)
        {
            eptr[i] = cells[i].first_node_idx;
            for (auto j = 0; j < cells[i].n_nodes; j++)
            {
                eind[eptr[i] + j] = cell_node_conn[eptr[i] + j];
            }
        }
        eptr[n_cells] = n_cell_node_idx;

        // Execute METIS partitioning routine
        // Check return value for errors
        int metis_obj_val;
        int metis_return_val = METIS_PartMeshNodal((idx_t *)&n_cells,
                                                   (idx_t *)&n_nodes,
                                                   eptr.data(),
                                                   eind.data(),
                                                   nullptr,
                                                   nullptr,
                                                   &n_parts,
                                                   nullptr,
                                                   nullptr,
                                                   &metis_obj_val,
                                                   cpart.data(),
                                                   npart.data());
        if (metis_return_val != 1)
        {
            std::string message = "METIS_PartMeshNodal() returned with:" + std::to_string(metis_return_val) + ", exiting\n";
            Logger::GetInstance().Error(message, __FILE__, __LINE__);
        }

        return std::make_pair(cpart, npart);
    }
#endif // SFEM_USE_METIS
}