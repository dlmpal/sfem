#include "gmsh.h"

namespace sfem::io
{
    //=============================================================================
    mesh::Mesh ReadGmsh(const std::string &filename)
    {
        std::ifstream file(filename);
        if (!file.is_open())
        {
            error::InvalidFileNameError(filename, __FILE__, __LINE__);
        }

        // Keep track of current file line
        int line_idx = 0;

        // Skip first 4 lines
        std::string buffer;
        for (auto i = 0; i < 6; i++)
        {
            file >> buffer;
            line_idx++;
        }

        // Read regions ($PhysicalNames)
        int n_regions;
        std::vector<mesh::Region> regions;
        file >> n_regions;
        line_idx++;
        for (auto i = 0; i < n_regions; i++)
        {
            std::string name;
            int dim, tag;
            file >> dim;
            file >> tag;
            file >> name;
            name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
            regions.push_back(mesh::Region(name, dim, tag));
            line_idx += 3;
        }

        // Skip two lines
        for (auto i = 0; i < 2; i++)
        {
            file >> buffer;
            line_idx++;
        }

        // Read nodal positions ($Nodes)
        int n_nodes;
        file >> n_nodes;
        line_idx++;
        std::vector<Float> xpts(n_nodes * 3);
        for (auto i = 0; i < n_nodes; i++)
        {
            file >> buffer;
            for (auto j = 0; j < 3; j++)
            {
                file >> xpts[i * 3 + j];
            }
            line_idx += 4;
        }

        // Skip two lines
        for (auto i = 0; i < 2; i++)
        {
            file >> buffer;
            line_idx++;
        }

        // Read elements ($Elements), but skip the connectivity
        int n_cells, first_node_idx = 0;
        file >> n_cells;
        line_idx++;
        std::vector<mesh::Cell> cells(n_cells);
        for (auto i = 0; i < n_cells; i++)
        {
            int id, gmsh_type, dim, physical_tag, elementary_tag;
            file >> id;
            file >> gmsh_type;
            file >> dim;
            file >> physical_tag;
            file >> elementary_tag;
            auto [type, order] = gmsh::GmshTypeToNative(id, gmsh_type);
            cells[i] = mesh::Cell(id - 1, type, order, physical_tag, first_node_idx);
            for (auto j = 0; j < cells[i].n_nodes; j++)
            {
                file >> buffer;
            }
            first_node_idx += cells[i].n_nodes;
        }

        // Close the file and re-open it
        // at the $Elements section
        file.close();
        file.open(filename);
        for (auto i = 0; i < line_idx; i++)
        {
            file >> buffer;
        }

        // Read the cell-to-node connectivity
        std::vector<int> cell_node_conn(first_node_idx);
        for (auto i = 0; i < n_cells; i++)
        {
            for (auto j = 0; j < 5; j++)
            {
                file >> buffer;
            }
            for (auto j = 0; j < cells[i].n_nodes; j++)
            {
                file >> cell_node_conn[cells[i].first_node_idx + j];
            }
        }

        // Gmsh numbering starts at 1
        for (auto i = 0; i < cell_node_conn.size(); i++)
        {
            cell_node_conn[i] -= 1;
        }

        return mesh::Mesh(n_cells, cells, cell_node_conn, n_nodes, 0, n_nodes, xpts, regions);
    }
}