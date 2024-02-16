#pragma once

#include "../common/config.h"

namespace sfem::mesh
{
    /// @brief Available cell types
    enum class CellType : int
    {
        INVALID = -1,
        POINT = 0,
        LINE = 1,
        TRIANGLE = 2,
        QUAD = 3,
        TET = 4,
        HEX = 5,
        PRISM = 6
    };

    /// @brief Number of nodes per cell type and order
    inline int GetCellNumNodes(CellType type, int order)
    {
        static const std::map<std::pair<CellType, int>, int> n_nodes_per_type = {
            {{CellType::POINT, 1}, 1},

            {{CellType::LINE, 1}, 2},
            {{CellType::LINE, 2}, 3},

            {{CellType::TRIANGLE, 1}, 3},
            {{CellType::TRIANGLE, 2}, 6},

            {{CellType::QUAD, 1}, 4},
            {{CellType::QUAD, 2}, 8},

            {{CellType::TET, 1}, 4},
            {{CellType::TET, 2}, 10},

            {{CellType::HEX, 1}, 8},
            {{CellType::HEX, 2}, 20},

            {{CellType::PRISM, 1}, 6}};

        if (n_nodes_per_type.count({type, order}) == 0)
        {
            return 0;
        }
        else
        {
            return n_nodes_per_type.at(std::make_pair(type, order));
        }
    }

    /// @brief POD t
    struct Cell
    {
        // Default constructor
        Cell() = default;

        /// @brief Create a Cell
        Cell(int id, CellType type, int order, int region_tag, int first_node_idx);

        /// @brief Global index
        int id = -1;

        /// @brief Type
        CellType type = CellType::INVALID;

        /// @brief Order
        int order = -1;

        /// @brief Integer tag of the Region to which the Cell belongs
        int region_tag = -1;

        /// @brief Number of nodes
        int n_nodes = -1;

        /// @brief The starting index for the Cell inside the the cell-to-node connecivity vector
        int first_node_idx = -1;
    };
}