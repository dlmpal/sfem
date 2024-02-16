#pragma once

#include "config.h"

namespace sfem
{
    struct Logger
    {
        /* Standard info message */
        static int const INFO = 0;

        /* Program will not abort,
        but normal execution is likely disturbed */
        static int const WARN = 1;

        /* Program will abort */
        static int const ERROR = -1;

        /* Process rank; */
        int proc_rank;

        /* Number of processes */
        int n_procs;

        static Logger &GetInstance(const std::string &appliaction_name = "SFEM_PROGRAM", int proc_rank = SFEM_ROOT, int n_procs = 1);

        void LogMessage(const std::string &message, int level = INFO, const std::string &file = "", int line = -1) const;

    private:
        Logger(const std::string &appliaction_name = "SFEM_PROGRAM", int proc_rank = SFEM_ROOT, int n_procs = 1);

        Logger(const Logger &) = delete;

        Logger &operator=(const Logger &) = delete;

        std::string LevelToString(int level) const;

        void Abort() const;

        static Logger *_logger;
    };
}