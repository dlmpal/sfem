#include "region.h"

namespace sfem::mesh
{
    //=============================================================================
    Region::Region(const std::string &name, int dim, int tag)
        : name(name), dim(dim), tag(tag)
    {
    }

}