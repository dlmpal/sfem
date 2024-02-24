#pragma once

#include "../common/config.h"

namespace sfem::mesh
{
  /// @brief POD
  struct Region
  {
    // Default constructor
    Region() = default;

    /// @brief Create a region
    Region(const std::string &name, int dim, int tag);

    /// @brief Name
    std::string name = {};

    /// @brief Physical dimension
    int dim = -1;

    /// @brief Numerical tag
    int tag = -1;
  };
}
