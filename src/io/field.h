#pragma once

#include "../core/field.h"

namespace sfem::io
{
    /// @brief Read field values from file
    void ReadFieldValues(const std::string &path, field::Field &field);

    /// @brief Write field values to file
    void WriteFieldValues(const std::string &path, const field::Field &field, bool assemble_global = true);
}