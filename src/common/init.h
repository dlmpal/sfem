#pragma once

#include "config.h"

namespace sfem
{
    /// @brief Should be called at the start of each sfem application. Returns process rank and number of processes
    /// @returns The process rank and the number of processes
    void Initialize(int *argc, char ***argv, const std::string &application_name = "SFEM_PROGRAM", const std::string &log_path = "");

    /// @brief Should be called at the end of each sfem application
    void Finalize();
}