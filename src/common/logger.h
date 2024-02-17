#pragma once

#include "config.h"

namespace sfem
{
    class Logger
    {
    public:
        enum class LogLevel : int
        {
            INFO = 0,
            WARN = 1,
            ERROR = 2
        };

        Logger(const Logger &) = delete;
        Logger &operator=(const Logger &) = delete;

        /// @brief Get the (static) Logger instance
        /// @note All input parameters only matter for initialization
        /// @param appliaction_name The name of the SFEM application
        /// @param path The path for the log file
        /// @param proc_rank MPI process rank
        /// @param n_procs Number of MPI processes
        static Logger &GetInstance(const std::string &appliaction_name = "", const std::string &path = "", int proc_rank = SFEM_ROOT, int n_procs = 1);

        ///@brief Get the MPI process rank
        int GetProcRank() const;

        /// @brief Get the number of MPI processes
        int GetNumProcs() const;

        /// @brief Log a message
        /// @param message The messaged to be logged
        /// @param level The LogLevel
        /// @param file The file from which the message originates (e.g __FILE__)
        /// @param line The line in the file (e.g __LINE__)
        void LogMessage(const std::string &message, LogLevel level = LogLevel::INFO, const std::string &file = "", int line = -1);

        /// @brief Log an info message
        /// @note See LogMessage
        void Info(const std::string &message);

        /// @brief Log a warning
        /// @note See LogMessage
        void Warn(const std::string &message, const std::string &file = "", int line = -1);

        /// @brief Log an error
        /// @note Also aborts the program
        /// @note See LogMessage
        void Error(const std::string &message, const std::string &file = "", int line = -1);

        /// @brief Abort the application
        void Abort() const;

        /// @brief Get the corresponding string for a LogLevel
        std::string LevelToString(LogLevel level) const;

    private:
        /// @brief Constructor
        /// @note See GetInstance
        Logger(const std::string &appliaction_name = "SFEM_PROGRAM", const std::string &path = "", int proc_rank = SFEM_ROOT, int n_procs = 1);

        int proc_rank;
        int n_procs;
        std::ofstream log_file;
    };
}