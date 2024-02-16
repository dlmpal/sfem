#pragma once

#include "solver.h"
#include "../time_integrators/time_integrator.h"

namespace sfem::solver
{
    struct DynamicSolver : Solver
    {

        /* Starting time */
        Float t_start;

        /* End time */
        Float t_final;

        /* Time step size */
        Float dt;

        /* Current time context */
        std::pair<int, Float> time;

        /* Time-stepping */
        time_integrator::TimeIntegrator *integrator;

        DynamicSolver(assembly::MonoFieldAssembler *assembler, Float t_start, Float t_final, Float dt, int integrator_type = time_integrator::IMPLICIT_EULER);

        ~DynamicSolver();

        void Run();

        void _Run();

        void PostRun();

        /* Get number of time steps */
        int GetNumTimeSteps();

        /* Reset time context */
        void ResetTime();

        /* Get the current time context */
        std::pair<int, Float> GetCurrentTime();

        /* Increment time context */
        void IncrementTime();
    };
}