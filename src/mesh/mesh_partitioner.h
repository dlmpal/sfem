#pragma once

#include "mesh.h"

namespace sfem::mesh
{
    class MeshPartitioner
    {
    public:
        MeshPartitioner(const Mesh *mesh, int n_parts);
        void PartitionMesh();
        void WritePartition(const std::string &path);

    protected:
        virtual std::pair<std::vector<int>, std::vector<int>> CreatePartition() = 0;
        void AssignCells(std::vector<int> epart);
        void AssignNodes(std::vector<int> npart);
        void RenumberNodes();
        void WriteCellPartition(const std::string &path);
        void WriteNodePartition(const std::string &path);

        /// @brief Mesh to be partioned
        const Mesh *mesh;

        /// @brief Number of partitions
        int n_parts;

        /// @brief Cell partition
        std::vector<std::vector<int>> cells_per_part;
        std::vector<int> n_cell_node_idx_per_part;

        /// @brief Node partition
        std::vector<std::vector<int>> nodes_owned_per_part;
        std::vector<std::vector<int>> nodes_ghost_per_part;
        std::vector<std::vector<int>> nodes_owned_per_part_renumbered;
        std::vector<std::vector<int>> nodes_ghost_per_part_renumbered;
    };

    /// @brief Constructs a MeshPartioner given the type, e.g METIS
    MeshPartitioner *CreateMeshPartitioner(const std::string &type, const mesh::Mesh *mesh, int n_parts);

#ifdef SFEM_USE_METIS
    /// @brief Partition the Mesh using METIS
    class METISPartitioner : public MeshPartitioner
    {
    public:
        METISPartitioner(const Mesh *mesh, int n_parts);

    private:
        std::pair<std::vector<int>, std::vector<int>> CreatePartition() override;
    };
#endif
}