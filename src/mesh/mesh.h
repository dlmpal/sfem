#pragma once

#include "cell.h"
#include "region.h"

namespace sfem::mesh
{
    /// @brief Indexing types:
    enum class IndexingType : int
    {
        LOCAL = 0,
        GLOBAL = 1,
        RENUMBERED = 2
    };

    class Mesh
    {
    public:
        /// @brief Create a Mesh (serial)
        /// @note The required index mappings are also (trivially) generated
        /// @param cells Cell vector
        /// @param cell_node_conn Cell-to-node connectivity vector
        /// @param xpts Nodal positions
        /// @param regions Region vector
        Mesh(std::vector<Cell> cells, std::vector<int> cell_node_conn, std::vector<Float> xpts, std::vector<Region> regions);

        /// @brief Create a mesh (parallel/partioned)
        /// @param n_cells_global Global number of cells
        /// @param cells Cell vector
        /// @param cell_node_conn Cell-to-node connectivity vector
        /// @param n_nodes_owned Number of owned nodes for this process
        /// @param n_nodes_ghost Number of ghost nodes for this process
        /// @param n_nodes_global Global number of nodes
        /// @param xpts Nodal positions
        /// @param regions Region vector
        /// @param node_local_to_global Node local-to-global mapping
        /// @param node_global_to_local Node global-to-local mapping
        Mesh(int n_cells_global, std::vector<Cell> cells, std::vector<int> cell_node_conn,
             int n_nodes_owned, int n_nodes_ghost, int n_nodes_global, std::vector<Float> xpts, std::vector<Region> regions,
             std::pair<std::vector<int>, std::vector<int>> node_local_to_global = {},
             std::pair<std::unordered_map<int, int>, std::unordered_map<int, int>> node_global_to_local = {}); // Distributed

        /// @brief Print information about the Mesh
        void Info() const;

        /// @brief Get the physical dimension
        int GetDim() const;

        /// @brief Get the local number of cells
        int GetNumCells() const;

        /// @brief Get the global number of cells
        int GetNumCellsGlobal() const;

        /// @brief Get the size of the cell-to-node connectivity vector
        int GetSizeCellNodeConn() const;

        /// @brief Get the local number of nodes
        /// @note The ghost nodes are also included
        int GetNumNodes() const;

        /// @brief Get the number owned nodes for this process
        int GetNumNodesOwned() const;

        /// @brief Get the number of ghost nodes for this process
        int GetNumNodesGhost() const;

        /// @brief Get the global number of nodes
        int GetNumNodesGlobal() const;

        /// @brief Get the number of regions
        int GetNumRegions() const;

        /// @brief Get a reference to the Cell vector
        const std::vector<Cell> &GetCells() const;

        /// @brief Get a reference to the cell-to-node connectivity vector
        const std::vector<int> &GetCellNodeConn() const;

        /// @brief Get a reference to the nodal positions vector
        const std::vector<Float> &GetXpts() const;

        /// @brief Get a reference to the Region vector
        const std::vector<Region> &GetRegions() const;

        /// @brief Get the desired local-to-global mapping for the nodes
        const std::vector<int> &GetNodeLocalToGlobalMapping(IndexingType indexing) const;

        /// @brief Get the desired global-to-local mapping for the nodes
        const std::unordered_map<int, int> &GetNodeGlobalToLocalMapping(IndexingType indexing) const;

        /// @brief Map nodes from local to desired indexing
        /// @param n_nodes Number of nodes
        /// @param indexing Desired indexing
        /// @param nodes Nodes
        void MapLocalNodes(int n_nodes, IndexingType indexing, int nodes[]) const;

        /// @brief Returns true if a node is owned locally
        /// @note Only works for local indexing
        bool IsNodeOwned(int node) const;

        /// @brief Get the owned nodes for this process in the desired indexing
        std::vector<int> GetOwnedNodes(IndexingType indexing) const;

        /// @brief Get the ghost nodes for this process in the desired indexing
        std::vector<int> GetGhostNodes(IndexingType indexing) const;

        /// @brief Get a Region by its name
        Region GetRegionByName(const std::string &region_name) const;

        /// @brief Get the Cells for a Region, given its name
        std::vector<Cell> GetRegionCells(const std::string &region_name) const;

        /// @brief Get the nodes for a Region, given its name, in the desired indexing
        std::vector<int> GetRegionNodes(const std::string &region_name, IndexingType indexing) const;

        /// @brief Get the nodes and their positions for a certain cell, in the desired indexing
        std::pair<std::array<int, SFEM_MAX_CELL_NODES>, std::array<Float, SFEM_MAX_CELL_NODES * 3>>
        GetCellNodes(const Cell &cell, IndexingType indexing = IndexingType::LOCAL) const;

    private:
        /// @brief Physical dimension
        int dim;

        /// @brief Global number of Cells
        int n_cells_global;

        /// @brief Vector containing the local Cells
        std::vector<Cell> cells;

        /// @brief Vector describing the cell-to-node connectivity
        /// @note For a Cell with local index i, its corresponding nodes are
        /// found by accessing cell_node_conn[i + k], where k goes from 0 to Cell.n_nodes
        std::vector<int> cell_node_conn;

        /// @brief Number of nodes owned by this process
        int n_nodes_owned;

        /// @brief Number of ghost nodes for this process
        int n_nodes_ghost;

        /// @brief Global number of nodes
        int n_nodes_global;

        /// @brief Nodal positions
        std::vector<Float> xpts;

        /// @brief Vector containing the Regions
        std::vector<Region> regions;

        /// @brief Map nodes from local to global (first) or renumbered (second) indexing
        std::pair<std::vector<int>, std::vector<int>> node_local_to_global;

        /// @brief Map nodes from global (first) /renumbered (second) to local indexing
        std::pair<std::unordered_map<int, int>, std::unordered_map<int, int>> node_global_to_local;
    };

}