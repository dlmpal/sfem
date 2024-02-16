#pragma once

#include "../common/config.h"

namespace sfem::mesh
{
  struct Region
  {
    Region() = default;
    Region(const std::string &name, int dim, int tag);

    // Data
    std::string name = {};
    int dim = -1;
    int tag = -1;
  };
}
