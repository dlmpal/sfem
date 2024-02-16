#include "mesh.h"
#include "../common/error.h"

namespace sfem::io
{
    //=============================================================================
    std::tuple<int, int, std::pair<std::vector<int>, std::vector<int>>, std::pair<std::unordered_map<int, int>, std::unordered_map<int, int>>>
    ReadNodePartition(const std::string &path)
    {
        int n_procs = Logger::GetInstance().n_procs;
        int proc_rank = Logger::GetInstance().proc_rank;

        std::ifstream file(path + "/NodePartition");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/NodePartition", __FILE__, __LINE__);
        }

        // Check the number of partitions is equal to the number of processes
        int n_parts;
        file >> n_parts;
        if (n_parts != n_procs)
        {
            error::InvalidSizeError(n_parts, Logger::GetInstance().n_procs, __FILE__, __LINE__);
        }

        // Read the number of owned+ghost nodes per process
        int n_nodes_owned_per_part[n_parts];
        int n_nodes_ghost_per_part[n_parts];
        for (auto i = 0; i < n_parts; i++)
        {
            file >> n_nodes_owned_per_part[i];
            file >> n_nodes_ghost_per_part[i];
        }
        int n_nodes_owned = n_nodes_owned_per_part[proc_rank];
        int n_nodes_ghost = n_nodes_ghost_per_part[proc_rank];
        int n_nodes = n_nodes_owned + n_nodes_ghost;

        // Skip to this process' section of the file
        int offset = 0;
        for (auto i = 0; i < proc_rank; i++)
        {
            offset += n_nodes_owned_per_part[i] + n_nodes_ghost_per_part[i];
        }
        int buffer;
        for (auto i = 0; i < offset * 2; i++)
        {
            file >> buffer;
        }

        // Read the index mappings
        std::vector<int> local_to_global(n_nodes);
        std::vector<int> local_to_renumbered(n_nodes);
        int global, renumbered;
        for (auto i = 0; i < n_nodes; i++)
        {
            file >> global;
            file >> renumbered;
            local_to_global[i] = global;
            local_to_renumbered[i] = renumbered;
        }

        std::unordered_map<int, int> global_to_local;
        std::unordered_map<int, int> renumbered_to_local;
        for (auto i = 0; i < n_nodes; i++)
        {
            int global = local_to_global.at(i);
            int renumbered = local_to_renumbered.at(i);
            global_to_local[global] = i;
            renumbered_to_local[renumbered] = i;
        }

        return std::make_tuple(n_nodes_owned, n_nodes_ghost, std::make_pair(local_to_global, local_to_renumbered), std::make_pair(global_to_local, renumbered_to_local));
    }
    //=============================================================================
    std::tuple<int, int, std::vector<int>> ReadCellPartition(const std::string &path)
    {
        int n_procs = Logger::GetInstance().n_procs;
        int proc_rank = Logger::GetInstance().proc_rank;

        std::ifstream file(path + "/CellPartition");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/CellPartition", __FILE__, __LINE__);
        }

        // Check the number of partitions is equal to the number of processes
        int n_parts;
        file >> n_parts;
        if (n_parts != n_procs)
        {
            error::InvalidSizeError(n_parts, Logger::GetInstance().n_procs, __FILE__, __LINE__);
        }

        // Read the number of cells per process
        int n_cells_per_part[n_parts], size_cell_node_conn_per_part[n_parts];
        for (auto i = 0; i < n_parts; i++)
        {
            file >> n_cells_per_part[i];
            file >> size_cell_node_conn_per_part[i];
        }
        int n_cells = n_cells_per_part[proc_rank];
        int size_cell_node_conn = size_cell_node_conn_per_part[proc_rank];

        // Skip to this process' section of the file
        int offset = 0;
        for (auto i = 0; i < proc_rank; i++)
        {
            offset += n_cells_per_part[i];
        }
        int buffer;
        for (auto i = 0; i < offset; i++)
        {
            file >> buffer;
        }

        // Read the index mapping
        std::vector<int> local_to_global(n_cells);
        for (auto i = 0; i < n_cells; i++)
        {
            file >> local_to_global[i];
        }

        // Sort so that the global cell indices are in ascending order
        std::sort(local_to_global.begin(), local_to_global.end());

        return std::make_tuple(n_cells, size_cell_node_conn, local_to_global);
    }
    //=============================================================================
    std::tuple<int, std::vector<Float>> ReadNodes(const std::string &path, const std::unordered_map<int, int> &global_to_local = {})
    {
        int n_procs = Logger::GetInstance().n_procs;

        std::ifstream file(path + "/Nodes");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Nodes", __FILE__, __LINE__);
        }

        // Read the global number of nodes
        int n_nodes_global;
        file >> n_nodes_global;

        // Number of owned + ghost nodes
        int n_nodes;
        if (n_procs == 1)
        {
            n_nodes = n_nodes_global;
        }
        else
        {
            n_nodes = global_to_local.size();
        }

        // Read nodal positions
        std::vector<Float> xpts(n_nodes * 3);
        Float xpt[3];
        for (auto i = 0; i < n_nodes_global; i++)
        {
            for (auto j = 0; j < 3; j++)
            {
                file >> xpt[j];
            }
            if (n_procs == 1 || global_to_local.count(i) == 1)
            {
                int local = n_procs == 1 ? i : global_to_local.at(i);
                for (auto j = 0; j < 3; j++)
                {
                    xpts[local * 3 + j] = xpt[j];
                }
            }
        }

        return std::make_tuple(n_nodes_global, xpts);
    }
    //=============================================================================
    std::tuple<int, std::vector<mesh::Cell>, std::vector<int>>
    ReadCells(const std::string &path, int _size_cell_node_conn = 0, const std::vector<int> cell_local_to_global = {}, const std::unordered_map<int, int> &node_global_to_local = {})
    {
        int n_procs = Logger::GetInstance().n_procs;

        std::ifstream file(path + "/Cells");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Cells", __FILE__, __LINE__);
        }

        // Read the global number of cells and the global cell-to-node connectivity size
        int n_cells_global, size_cell_node_conn_global;
        file >> n_cells_global;
        file >> size_cell_node_conn_global;

        // Number of owned cells and cell-to-node connectivity size
        int n_cells, size_cell_node_conn;
        if (n_procs == 1)
        {
            n_cells = n_cells_global;
            size_cell_node_conn = size_cell_node_conn_global;
        }
        else
        {
            n_cells = cell_local_to_global.size();
            size_cell_node_conn = _size_cell_node_conn;
        }

        // Read cells and cell to node connectivity
        std::vector<mesh::Cell> cells(n_cells);
        std::vector<int> cell_node_conn(size_cell_node_conn);
        int first_node_idx = 0;
        int cell_idx = 0;
        for (auto i = 0; i < n_cells_global; i++)
        {
            int id, type, order, region_tag;
            file >> id;
            file >> type;
            file >> order;
            file >> region_tag;
            mesh::Cell cell(id, (mesh::CellType)type, order, region_tag, first_node_idx);

            // Read cell nodes
            int cell_nodes[cell.n_nodes];
            for (auto j = 0; j < cell.n_nodes; j++)
            {
                file >> cell_nodes[j];
            }

            // For distributed meshes, add only if cell is owned by this process
            if (n_procs == 1 || cell_local_to_global[cell_idx] == i)
            {
                cells[cell_idx] = cell;
                for (auto j = 0; j < cell.n_nodes; j++)
                {
                    cell_node_conn[first_node_idx] = cell_nodes[j];
                    first_node_idx++;
                }
                cell_idx++;
            }
        }

        // For distributed meshes, map the nodes to local indexing
        if (n_procs > 1)
        {
            for (auto i = 0; i < cell_node_conn.size(); i++)
            {
                cell_node_conn[i] = node_global_to_local.at(cell_node_conn[i]);
            }
        }

        return std::make_tuple(n_cells_global, cells, cell_node_conn);
    }
    //=============================================================================
    std::vector<mesh::Region> ReadRegions(const std::string &path)
    {
        std::ifstream file(path + "/Regions");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Regions", __FILE__, __LINE__);
        }

        // Read the number of regions
        int n_regions;
        file >> n_regions;

        // Read regions
        std::vector<mesh::Region> regions;
        for (auto i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> name;
            file >> dim;
            file >> tag;
            regions.push_back(mesh::Region(name, dim, tag));
        }

        return regions;
    }
    //=============================================================================
    mesh::Mesh ReadMesh(const std::string &path)
    {
        if (Logger::GetInstance().n_procs == 1)
        {
            auto [n_nodes, xpts] = ReadNodes(path);
            auto [n_cells, cells, cell_node_conn] = ReadCells(path);
            auto regions = ReadRegions(path);
            return mesh::Mesh(cells, cell_node_conn, xpts, regions);
        }
        else
        {
            auto [n_nodes_owned, n_nodes_ghost, node_local_to_global, node_global_to_local] = ReadNodePartition(path);
            auto [n_cells, size_cell_node_conn, cell_local_to_global] = ReadCellPartition(path);
            auto [n_nodes_global, xpts] = ReadNodes(path, node_global_to_local.first);
            auto [n_cells_global, cells, cell_node_conn] = ReadCells(path, size_cell_node_conn, cell_local_to_global, node_global_to_local.first);
            auto regions = ReadRegions(path);
            return mesh::Mesh(n_cells_global, cells, cell_node_conn,
                              n_nodes_owned, n_nodes_ghost, n_nodes_global, xpts, regions,
                              node_local_to_global, node_global_to_local);
        }
    }
    //=============================================================================
    void WriteNodes(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/Nodes");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Nodes", __FILE__, __LINE__);
        }
        file << mesh.GetNumNodes() << "\n";
        const auto &xpts = mesh.GetXpts();
        for (auto i = 0; i < mesh.GetNumNodes(); i++)
        {
            file << xpts[i * 3 + 0] << " ";
            file << xpts[i * 3 + 1] << " ";
            file << xpts[i * 3 + 2] << "\n";
        }
    }
    //=============================================================================
    void WriteCells(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/Cells");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Cells", __FILE__, __LINE__);
        }
        file << mesh.GetNumCells() << "\n";
        file << mesh.GetSizeCellNodeConn() << "\n";
        for (const auto &cell : mesh.GetCells())
        {
            file << cell.id << " " << (int)cell.type << " " << cell.order << " " << cell.region_tag << " ";
            auto [nodes, _] = mesh.GetCellNodes(cell);
            for (auto j = 0; j < cell.n_nodes; j++)
            {
                file << nodes[j] << " ";
            }
            file << "\n";
        }
    }
    //=============================================================================
    void WriteRegions(const std::string &path, const mesh::Mesh &mesh)
    {
        std::ofstream file(path + "/Regions");
        if (!file.is_open())
        {
            error::InvalidFileNameError(path + "/Regions", __FILE__, __LINE__);
        }
        file << mesh.GetNumRegions() << "\n";
        for (const auto &region : mesh.GetRegions())
        {
            file << region.name << " " << region.dim << " " << region.tag << "\n";
        }
    }
    //=============================================================================
    void WriteMesh(const std::string &path, const mesh::Mesh &mesh)
    {
        WriteNodes(path, mesh);
        WriteCells(path, mesh);
        WriteRegions(path, mesh);
    }
}