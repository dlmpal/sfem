#include "init.h"
#include "logger.h"
#include "../third_party/sfem_petsc.h"
#include "../third_party/sfem_slepc.h"

namespace sfem
{
    //=============================================================================
    void Initialize(int *argc, char ***argv, const std::string &application_name, const std::string &log_path)
    {
        int proc_rank = SFEM_ROOT, n_procs = 1;
        if (argc != nullptr)
        {
#ifdef SFEM_USE_PETSC
            PetscInitialize(argc, argv, nullptr, nullptr);
#endif
#ifdef SFEM_USE_MPI
            MPI_Comm_rank(SFEM_COMM_WORLD, &proc_rank);
            MPI_Comm_size(SFEM_COMM_WORLD, &n_procs);
#endif
#ifdef SFEM_USE_SLEPC
            SlepcInitialize(argc, argv, nullptr, nullptr);
#endif
        }
        Logger::GetInstance(application_name, log_path, proc_rank, n_procs);
    }
    //=============================================================================
    void Finalize()
    {
#ifdef SFEM_USE_PETSC
        PetscBool petsc_init;
        PetscInitialized(&petsc_init);
        if (petsc_init)
        {
            PetscFinalize();
        }
#endif
    }

}