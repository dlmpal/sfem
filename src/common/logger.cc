#include "logger.h"

namespace sfem
{
    //=============================================================================
    Logger::Logger(const std::string &appliaction_name, int proc_rank, int n_procs)
    {
        this->proc_rank = proc_rank;
        this->n_procs = n_procs;

        if (proc_rank == SFEM_ROOT)
        {
            std::string message = "Application: " + appliaction_name + "\n";
            message += "Number of processes launched: " + std::to_string(n_procs) + "\n";
            LogMessage(message);
        }
    }
    //=============================================================================
    Logger &Logger::GetInstance(const std::string &appliaction_name, int proc_rank, int n_procs)
    {
        static Logger *_logger = new Logger(appliaction_name, proc_rank, n_procs);
        return *_logger;
    }
    //=============================================================================
    void Logger::Abort() const
    {
#ifdef SFEM_USE_MPI
        MPI_Abort(SFEM_COMM_WORLD, ERROR);
#else
        exit(-1);
#endif
    }
    //=============================================================================
    std::string Logger::LevelToString(int level) const
    {
        switch (level)
        {
        case INFO:
            return std::string("INFO");

        case WARN:
            return std::string("WARNING");

        default:
            return std::string("ERROR");
        }
    }
    //=============================================================================
    void Logger::LogMessage(const std::string &message, int level, const std::string &file, int line) const
    {
        std::string _message = "#=============================================================================#\n";
        _message += "SFEM-PROCESS-" + std::to_string(proc_rank) + "-" + LevelToString(level) + "\n";
        if (level == ERROR)
            _message += "File: " + file + "\nLine: " + std::to_string(line) + "\n";
        _message += message;
        _message += "#=============================================================================#\n";

        std::cout << _message;
        if (level == ERROR)
            Abort();
    }
}