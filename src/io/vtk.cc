#include "vtk.h"

namespace sfem::io
{
    //=============================================================================
    void WriteVTK(const std::string &filename, const mesh::Mesh &mesh, const std::vector<field::Field> &fields)
    {
        std::ofstream file(filename);
        if (!file.is_open())
        {
            error::InvalidFileNameError(filename, __FILE__, __LINE__);
        }

        // File header
        file << "# vtk DataFile Version 2.0\n";
        file << "SFEM\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Points
        file << "POINTS " << mesh.GetNumNodes() << " float\n";
        const auto &xpts = mesh.GetXpts();
        for (auto i = 0; i < mesh.GetNumNodes(); i++)
        {
            file << xpts[i * 3 + 0] << " " << xpts[i * 3 + 1] << " " << xpts[i * 3 + 2] << "\n";
        }

        // Cells
        int vtk_size = mesh.GetSizeCellNodeConn() + mesh.GetNumCells();
        file << "CELLS " << mesh.GetNumCells() << " " << vtk_size << "\n";
        for (const auto &cell : mesh.GetCells())
        {
            file << cell.n_nodes << " ";
            auto [nodes, _] = mesh.GetCellNodes(cell);
            vtk::CellNodeOrderingVTK(cell, nodes.data());
            for (auto j = 0; j < cell.n_nodes; j++)
            {
                file << nodes[j] << " ";
            }
            file << "\n";
        }

        // Cell types
        file << "CELL_TYPES " << mesh.GetNumCells() << "\n";
        for (const auto &cell : mesh.GetCells())
        {
            file << vtk::CellTypeVTK(cell) << "\n";
        }

        // Field values
        if (fields.size() > 0)
        {
            file << "POINT_DATA " << mesh.GetNumNodes() << "\n";
            for (const auto &field : fields)
            {
                int n_vars = field.GetNumVars();
                if (n_vars > 1)
                {
                    file << "VECTORS "
                         << " " << field.GetName() << " float\n";
                }
                if (n_vars == 1)
                {
                    file << "SCALARS "
                         << " " << field.GetName() << " float\n";
                    file << "LOOKUP_TABLE default\n";
                }
                const auto &values = field.GetValues();
                for (auto i = 0; i < mesh.GetNumNodes(); i++)
                {
                    for (auto j = 0; j < n_vars; j++)
                    {
                        file << values[i * n_vars + j] << " ";
                    }
                    if (n_vars == 2)
                    {
                        file << "0.0 ";
                    }
                    file << "\n";
                }
            }
        }
    }
}