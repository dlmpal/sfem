#pragma once

#include "config.h"
#include "logger.h"

namespace sfem
{
    /// @brief Should be called at the start of each sfem application. Returns process rank and number of processes
    /// @returns The process rank and the number of processes
    std::pair<int, int> Initialize(std::string application_name, int *argc, char ***argv);

    /// @brief Should be called at the end of each sfem application
    void Finalize();
}