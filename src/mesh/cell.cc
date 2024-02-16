#include "cell.h"
#include "../common/error.h"

namespace sfem::mesh
{
    //=============================================================================
    Cell::Cell(int id, CellType type, int order, int region_tag, int first_node_idx)
        : id(id), type(type), order(order), region_tag(region_tag), first_node_idx(first_node_idx)
    {
        n_nodes = GetCellNumNodes(type, order);
        if (n_nodes == 0)
        {
            error::InvalidCellError(id, (int)type, order, __FILE__, __LINE__);
        }
    }
}