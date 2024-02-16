#pragma once

#include "../mesh/mesh.h"

namespace sfem::la
{
    /// @brief
    /// @param mesh
    /// @param n_vars
    /// @return
    std::pair<std::vector<int>, std::vector<int>> SparsityPattern(const mesh::Mesh *mesh, int n_vars);
}