#include "mesh.h"
#include "../common/logger.h"
#include "../common/error.h"

namespace sfem::mesh
{
    //=============================================================================
    Mesh::Mesh(std::vector<Cell> cells, std::vector<int> cell_node_conn,
               std::vector<Float> xpts, std::vector<Region> regions)
    {
        // Cells
        this->n_cells_global = cells.size();
        this->cells = cells;
        this->cell_node_conn = cell_node_conn;

        // Nodes
        this->n_nodes_owned = xpts.size() / 3;
        this->n_nodes_ghost = 0;
        this->n_nodes_global = this->n_nodes_owned;
        this->xpts = xpts;

        // Regions
        this->regions = regions;
        this->dim = 0;
        for (auto region : regions)
        {
            if (region.dim > dim)
            {
                dim = region.dim;
            }
        }

        // Create the index mappings
        // For serial meshes the mappings are trivial
        this->node_local_to_global.first.resize(this->n_nodes_owned);
        this->node_local_to_global.second.resize(this->n_nodes_owned);
        for (auto i = 0; i < this->n_nodes_owned; i++)
        {
            this->node_local_to_global.first[i] = i;
            this->node_local_to_global.second[i] = i;
            this->node_global_to_local.first.insert({i, i});
            this->node_global_to_local.second.insert({i, i});
        }
    }
    //=============================================================================
    Mesh::Mesh(int n_cells_global, std::vector<Cell> cells, std::vector<int> cell_node_conn,
               int n_nodes_owned, int n_nodes_ghost, int n_nodes_global, std::vector<Float> xpts, std::vector<Region> regions,
               std::pair<std::vector<int>, std::vector<int>> node_local_to_global,
               std::pair<std::unordered_map<int, int>, std::unordered_map<int, int>> node_global_to_local)
    {
        // Cells
        this->n_cells_global = n_cells_global;
        this->cells = cells;
        this->cell_node_conn = cell_node_conn;

        // Nodes
        this->n_nodes_global = n_nodes_global;
        this->n_nodes_owned = n_nodes_owned;
        this->n_nodes_ghost = n_nodes_ghost;
        this->xpts = xpts;
        if (GetNumNodes() != xpts.size() / 3)
        {
            error::InvalidSizeError(GetNumNodes(), xpts.size() / 3, __FILE__, __LINE__);
        }

        // Regions
        this->regions = regions;
        this->dim = 0;
        for (auto region : regions)
        {
            if (region.dim > this->dim)
            {
                this->dim = region.dim;
            }
        }

        // Index mappings
        this->node_local_to_global = node_local_to_global;
        this->node_global_to_local = node_global_to_local;
    }
    //=============================================================================
    void Mesh::Info() const
    {
        std::string message = "Mesh Info\n";
        message += "*Number of Nodes: " + std::to_string(GetNumNodes()) + "\n";
        message += "*Number of Cells: " + std::to_string(GetNumCells()) + "\n";
        message += "*Number of Regions: " + std::to_string(GetNumRegions()) + "\n";
        message += "*Regions: \n";
        for (auto region : regions)
            message += "\t|-" + region.name + "\n";
        Logger::GetInstance().LogMessage(message, Logger::INFO);
    }
    //=============================================================================
    int Mesh::GetDim() const
    {
        return dim;
    }
    //=============================================================================
    int Mesh::GetNumCells() const
    {
        return cells.size();
    }
    //=============================================================================
    int Mesh::GetNumCellsGlobal() const
    {
        return n_cells_global;
    }
    //=============================================================================
    int Mesh::GetSizeCellNodeConn() const
    {
        return cell_node_conn.size();
    }
    //=============================================================================
    int Mesh::GetNumNodes() const
    {
        return n_nodes_owned + n_nodes_ghost;
    }
    //=============================================================================
    int Mesh::GetNumNodesOwned() const
    {
        return n_nodes_owned;
    }
    //=============================================================================
    int Mesh::GetNumNodesGhost() const
    {
        return n_nodes_ghost;
    }
    //=============================================================================
    int Mesh::GetNumNodesGlobal() const
    {
        return n_nodes_global;
    }
    //=============================================================================
    int Mesh::GetNumRegions() const
    {
        return regions.size();
    }
    //=============================================================================
    const std::vector<Float> &Mesh::GetXpts() const
    {
        return xpts;
    }
    //=============================================================================
    const std::vector<Cell> &Mesh::GetCells() const
    {
        return cells;
    }
    //=============================================================================
    const std::vector<int> &Mesh::GetCellNodeConn() const
    {
        return cell_node_conn;
    }
    //=============================================================================
    const std::vector<Region> &Mesh::GetRegions() const
    {
        return regions;
    }
    //=============================================================================
    const std::vector<int> &Mesh::GetNodeLocalToGlobalMapping(IndexingType indexing) const
    {
        if (indexing == IndexingType::GLOBAL)
        {
            return node_local_to_global.first;
        }
        else if (indexing == IndexingType::RENUMBERED)
        {
            return node_local_to_global.second;
        }
        else
        {
            /// @todo: INVALID INDEXING ERROR;
        }
    }
    //=============================================================================
    const std::unordered_map<int, int> &Mesh::GetNodeGlobalToLocalMapping(IndexingType indexing) const
    {
        if (indexing == IndexingType::GLOBAL)
        {
            return node_global_to_local.first;
        }
        else if (indexing == IndexingType::RENUMBERED)
        {
            return node_global_to_local.second;
        }
        else
        {
            /// @todo: INVALID INDEXING ERROR;
        }
    }
    //=============================================================================
    void Mesh::MapLocalNodes(int n_nodes, IndexingType indexing, int nodes[]) const
    {
        const std::vector<int> *mapping = nullptr;
        if (indexing == IndexingType::LOCAL)
        {
            return;
        }
        else if (indexing == IndexingType::GLOBAL)
        {
            mapping = &node_local_to_global.first;
        }
        else if (indexing == IndexingType::RENUMBERED)
        {
            mapping = &node_local_to_global.second;
        }
        else
        {
            /// @todo: INVALID INDEXING ERROR;
        }
        for (auto i = 0; i < n_nodes; i++)
        {
            nodes[i] = (*mapping)[nodes[i]];
        }
    }
    //=============================================================================
    bool Mesh::IsNodeOwned(int node) const
    {
        if (node >= n_nodes_owned)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
    //=============================================================================
    std::vector<int> Mesh::GetOwnedNodes(IndexingType indexing) const
    {
        std::vector<int> owned_nodes(n_nodes_owned);
        for (auto i = 0; i < n_nodes_owned; i++)
        {
            owned_nodes[i] = i;
        }

        MapLocalNodes(owned_nodes.size(), indexing, owned_nodes.data());
        return owned_nodes;
    }
    //=============================================================================
    std::vector<int> Mesh::GetGhostNodes(IndexingType indexing) const
    {
        std::vector<int> ghost_nodes(n_nodes_ghost);
        for (auto i = 0; i < n_nodes_ghost; i++)
        {
            ghost_nodes[i] = i + n_nodes_owned;
        }
        MapLocalNodes(ghost_nodes.size(), indexing, ghost_nodes.data());
        return ghost_nodes;
    }
    //=============================================================================
    Region Mesh::GetRegionByName(const std::string &region_name) const
    {
        for (auto region : regions)
        {
            if (region.name == region_name)
            {
                return region;
            }
        }
        /// @todo: REGION ERROR
    }
    //=============================================================================
    std::vector<Cell> Mesh::GetRegionCells(const std::string &region_name) const
    {
        auto region = GetRegionByName(region_name);
        std::vector<Cell> region_cells;
        for (auto cell : cells)
        {
            if (cell.region_tag == region.tag)
            {
                region_cells.push_back(cell);
            }
        }
        return region_cells;
    }
    //=============================================================================
    std::vector<int> Mesh::GetRegionNodes(const std::string &region_name, IndexingType indexing) const
    {
        auto region_cells = GetRegionCells(region_name);
        std::vector<int> region_nodes;
        std::unordered_map<int, bool> is_node_included;
        for (auto cell : region_cells)
        {
            auto [nodes, _] = GetCellNodes(cell, indexing);
            for (auto i = 0; i < cell.n_nodes; i++)
            {
                if (is_node_included.count(nodes[i]) == 0)
                {
                    is_node_included[nodes[i]] = true;
                    region_nodes.push_back(nodes[i]);
                }
            }
        }
        return region_nodes;
    }
    //=============================================================================
    std::pair<std::array<int, SFEM_MAX_CELL_NODES>, std::array<Float, SFEM_MAX_CELL_NODES * 3>>
    Mesh::GetCellNodes(const Cell &cell, IndexingType indexing) const
    {
        std::array<int, SFEM_MAX_CELL_NODES> cell_nodes;
        std::array<Float, SFEM_MAX_CELL_NODES * 3> cell_xpts;
        for (auto i = 0; i < cell.n_nodes; i++)
        {
            cell_nodes[i] = cell_node_conn[cell.first_node_idx + i];
            for (auto j = 0; j < 3; j++)
            {
                cell_xpts[i * 3 + j] = xpts[cell_nodes[i] * 3 + j];
            }
        }
        MapLocalNodes(cell.n_nodes, indexing, cell_nodes.data());
        return std::make_pair(cell_nodes, cell_xpts);
    }
}