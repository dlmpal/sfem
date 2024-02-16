#pragma once

#include "logger.h"

/// @brief Commonly used errors
namespace sfem::error
{
    void InvalidFileNameError(const std::string &filename, const std::string &file, int line);
    void InvalidSizeError(int correct_size, int wrong_size, const std::string &file, int line);
    void InvalidCellError(int cell_id, int cell_type, int cell_order, const std::string &file, int line);
    void NegativeJacobianError(int cell_id, const std::string &file, int line);
    void InvalidFaceError(const std::string cell_type, int f_idx, const std::string &file, int line);
    void UnsupportedGmshTypeError(int gmsh_id, int gmsh_type, const std::string &file, int line);
}