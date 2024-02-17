#include "error.h"

namespace sfem::error
{
    //=============================================================================
    void InvalidFileNameError(const std::string &filename, const std::string &file, int line)
    {
        std::string message = "Error opening file: " + filename + "\n";
        Logger::GetInstance().Error(message, file, line);
    }
    //=============================================================================
    void InvalidSizeError(int correct_size, int wrong_size, const std::string &file, int line)
    {
        std::string message = "Expected size: " + std::to_string(correct_size) + "\n";
        message += "Got size: " + std::to_string(wrong_size) + "\n";
        Logger::GetInstance().Error(message, file, line);
    }
    //=============================================================================
    void InvalidCellError(int cell_id, int cell_type, int cell_order, const std::string &file, int line)
    {
        std::string message = "Cell " + std::to_string(cell_id) + " has invalid type (" + std::to_string(cell_type) + ") or order (" + std::to_string(cell_order) + ")\n";
        Logger::GetInstance().Error(message, file, line);
    }
    //=============================================================================
    void NegativeJacobianError(int cell_id, const std::string &file, int line)
    {
        std::string message = "Negative jacobian at cell: " + std::to_string(cell_id) + "\n";
        Logger::GetInstance().Error(message, file, line);
    }
    //=============================================================================
    void InvalidFaceError(std::string cell_type, int f_idx, const std::string &file, int line)
    {
        std::string message = "Invalid face index: " + std::to_string(f_idx) + "for " + cell_type + "\n";
        Logger::GetInstance().Error(message, file, line);
    }
    //=============================================================================
    void UnsupportedGmshTypeError(int gmsh_id, int gmsh_type, const std::string &file, int line)
    {
        std::string message = "Unsupported Gmsh cell type: " + std::to_string(gmsh_type) + "\n";
        Logger::GetInstance().Error(message, file, line);
    }

}