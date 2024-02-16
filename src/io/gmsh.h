#pragma once

#include "../common/error.h"
#include "../mesh/mesh.h"

namespace sfem::io
{
    /// @brief Read a Gmsh .msh2 file into sfem format
    /// @note This function is not designed for parallel execution
    mesh::Mesh ReadGmsh(const std::string &filename);
}

namespace sfem::io::gmsh
{
    /// @brief Get the sfem cell type and order from the corresponding Gmsh type
    inline std::pair<mesh::CellType, int> GmshTypeToNative(int gmsh_id, int gmsh_type)
    {
        static const std::map<int, std::pair<mesh::CellType, int>> from_gmsh = {
            {15, {mesh::CellType::POINT, 1}},

            {1, {mesh::CellType::LINE, 1}},
            {8, {mesh::CellType::LINE, 2}},

            {2, {mesh::CellType::TRIANGLE, 1}},
            {9, {mesh::CellType::TRIANGLE, 2}},

            {3, {mesh::CellType::QUAD, 1}},
            {16, {mesh::CellType::QUAD, 2}},

            {4, {mesh::CellType::TET, 1}},
            {11, {mesh::CellType::TET, 2}},

            {5, {mesh::CellType::HEX, 1}},
            {17, {mesh::CellType::HEX, 2}},

            {6, {mesh::CellType::PRISM, 1}},
        };
        if (from_gmsh.count(gmsh_type) == 0)
        {
            error::UnsupportedGmshTypeError(gmsh_id, gmsh_type, __FILE__, __LINE__);
        }
        return from_gmsh.at(gmsh_type);
    }
}