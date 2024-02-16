#pragma once

#include "../mesh/mesh.h"

namespace sfem::io
{
    /// @brief Read a Mesh from file
    /// @note For parallel execution, the NodePartition and CellPartition files are also required
    mesh::Mesh ReadMesh(const std::string &path);

    /// @brief Write a Mesh to file
    /// @note Not designed for parallel execution
    void WriteMesh(const std::string &path, const mesh::Mesh &mesh);
}