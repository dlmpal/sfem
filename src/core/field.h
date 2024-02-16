#pragma once

#include "../mesh/mesh.h"

/// @brief Field manipulation
namespace sfem::field
{
    class Field
    {
    public:
        /// @brief Create a Field
        /// @param name Name
        /// @param n_vars Number of variables
        /// @param mesh Mesh on which the field "lives"
        Field(const std::string &name, int n_vars, mesh::Mesh *mesh);

        /// @brief Get the name
        std::string GetName() const;

        /// @brief Get the number of variables
        int GetNumVars() const;

        /// @brief Get a pointer to the Mesh
        mesh::Mesh *GetMesh() const;

        /// @brief Get the local number of DOF
        /// @note Includes owned and ghost DOF
        int GetNumDof() const;

        /// @brief Get the number of owned DOF for this process
        int GetNumOwnedDof() const;

        /// @brief Get the number of ghost DOF for this process
        int GetNumGhostDof() const;

        /// @brief Get the global number of DOF
        int GetNumGlobalDof() const;

        /// @brief Get the owned DOF for this process, in the desired indexing
        std::vector<int> GetOwnedDof(mesh::IndexingType indexing = mesh::IndexingType::RENUMBERED) const;

        /// @brief Get the ghsot DOF for this process, in the desired indexing
        std::vector<int> GetGhostDof(mesh::IndexingType indexing = mesh::IndexingType::RENUMBERED) const;

        /// @brief Get the DOF corresponding to each given node
        /// @param[in] n_nodes Number of nodes
        /// @param[in] nodes Nodes
        /// @param[out] dof The corresponding DOF
        void GetNodeDof(int n_nodes, const int nodes[], int dof[]) const;

        /// @brief Get the DOF belonging to cell
        std::array<int, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD>
        GetCellDof(const mesh::Cell &cell, mesh::IndexingType indexing = mesh::IndexingType::RENUMBERED) const;

        /// @brief Get the values belonging to cell
        std::array<Float, SFEM_MAX_CELL_NODES * SFEM_MAX_N_VARS_FIELD>
        GetCellValues(const mesh::Cell &cell) const;

        /// @brief Assign fixed values to desired DOF
        void AddFixedDof(const std::string &region_name, int var, Float value);

        /// @brief Get the fixed DOF and their corresponding values
        std::pair<std::vector<int>, std::vector<Float>> GetFixedDof() const;

        /// @brief Clear all existing fixed dof
        void ClearFixedDof();

        /// @brief Set local DOF values
        void SetValues(const std::vector<Float> &values);

        /// @brief Get local DOF values
        const std::vector<Float> &GetValues() const;

    private:
        /// @brief Field name
        std::string name;

        /// @brief Number of Field variables
        int n_vars;

        /// @brief Mesh on which the Field "lives"
        mesh::Mesh *mesh;

        /// @brief Number of local DOF
        int n_dof;

        /// @brief Number of owned DOF
        int n_dof_owned;

        /// @brief Number of ghost DOF
        int n_dof_ghost;

        /// @brief Number of global DOF
        int n_dof_global;

        /// @brief Constrained DOF
        std::unordered_map<int, Float> fixed_dof;

        /// @brief Values corresponding to the local DOF
        std::vector<Float> values;
    };

    /// @brief Assemble all Field values to the root process
    std::vector<Float> AssembleGlobalValues(const Field &field);
}