#pragma once

#include "../common/error.h"
#include "../mesh/mesh.h"
#include "../core/field.h"

namespace sfem::io
{
    /// @brief Creates .vtk file for the given Mesh and Fields
    /// @note This function will create a file for each process that runs it.
    /// Therefore for each partition a separate file is created.
    void WriteVTK(const std::string &filename, const mesh::Mesh &mesh, const std::vector<field::Field> &fields);
}

namespace sfem::io::vtk
{
    /// @brief Get the corresponding VTK type for a sfem cell type
    inline int CellTypeVTK(const mesh::Cell &cell)
    {
        static const std::map<std::pair<mesh::CellType, int>, int> to_vtk = {
            {{mesh::CellType::POINT, 1}, 1},

            {{mesh::CellType::LINE, 1}, 3},
            {{mesh::CellType::LINE, 2}, 21},

            {{mesh::CellType::TRIANGLE, 1}, 5},
            {{mesh::CellType::TRIANGLE, 2}, 22},

            {{mesh::CellType::QUAD, 1}, 9},
            {{mesh::CellType::QUAD, 2}, 23},

            {{mesh::CellType::TET, 1}, 10},
            {{mesh::CellType::TET, 2}, 24},

            {{mesh::CellType::HEX, 1}, 12},
            {{mesh::CellType::HEX, 2}, 25},

            {{mesh::CellType::PRISM, 1}, 13}

        };
        if (to_vtk.count({cell.type, cell.order}) == 0)
        {
            error::InvalidCellError(cell.id, (int)cell.type, cell.order, __FILE__, __LINE__);
        }
        return to_vtk.at({cell.type, cell.order});
    }

    /// @brief Get the corresponding node ordering for VTK cells.
    /// @note For most sfem cell types no re-ordering is needed.
    inline void CellNodeOrderingVTK(const mesh::Cell &cell, int nodes[])
    {
        if (cell.type == mesh::CellType::TET and cell.order == 2)
        {
            int temp = nodes[9];
            nodes[9] = nodes[8];
            nodes[8] = temp;
        }
    }
}