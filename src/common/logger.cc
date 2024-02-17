#include "logger.h"

namespace sfem
{
    //=============================================================================
    Logger::Logger(const std::string &appliaction_name, const std::string &path, int proc_rank, int n_procs)
        : proc_rank(proc_rank), n_procs(n_procs)
    {
        if (!path.empty())
        {
            log_file = std::ofstream(path + "_" + std::to_string(proc_rank));
        }

        std::string message = "Application: " + appliaction_name + "\n";
        message += "Number of processes launched: " + std::to_string(n_procs) + "\n";
        LogMessage(message);
    }
    //=============================================================================
    Logger &Logger::GetInstance(const std::string &appliaction_name, const std::string &path, int proc_rank, int n_procs)
    {
        static Logger logger(appliaction_name, path, proc_rank, n_procs);
        return logger;
    }
    //=============================================================================
    int Logger::GetProcRank() const
    {
        return proc_rank;
    }
    //=============================================================================
    int Logger::GetNumProcs() const
    {
        return n_procs;
    }
    //=============================================================================
    void Logger::LogMessage(const std::string &message, LogLevel level, const std::string &file, int line)
    {
        std::string _message;
        _message += "#=============================================================================#\n";
        _message += "SFEM-PROCESS-" + std::to_string(proc_rank) + "-" + LevelToString(level) + "\n";
        if (level > LogLevel::INFO)
        {
            _message += "File: " + file + "\nLine: " + std::to_string(line) + "\n";
        }
        _message += message;
        _message += "#=============================================================================#\n";

        std::cout << _message;

        if (log_file.is_open())
        {
            log_file << _message;
        }

        if (level == LogLevel::ERROR)
        {
            Abort();
        }
    }
    //=============================================================================
    void Logger::Info(const std::string &message)
    {
        LogMessage(message);
    }
    //=============================================================================
    void Logger::Warn(const std::string &message, const std::string &file, int line)
    {
        LogMessage(message, LogLevel::WARN, file, line);
    }
    //=============================================================================
    void Logger::Error(const std::string &message, const std::string &file, int line)
    {
        LogMessage(message, LogLevel::ERROR, file, line);
    }
    //=============================================================================
    void Logger::Abort() const
    {
#ifdef SFEM_USE_MPI
        MPI_Abort(SFEM_COMM_WORLD, (int)LogLevel::ERROR);
#else
        exit(-1);
#endif
    }
    //=============================================================================
    std::string Logger::LevelToString(LogLevel level) const
    {
        switch (level)
        {
        case LogLevel::INFO:
            return std::string("INFO");
        case LogLevel::WARN:
            return std::string("WARNING");
        default:
            return std::string("ERROR");
        }
    }
}