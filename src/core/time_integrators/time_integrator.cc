#include "time_integrator.h"

namespace sfem::time_integrator
{
    TimeIntegrator::TimeIntegrator(assembly::MonoFieldAssembler *assembler, Float dt)
    {
        this->assembler = assembler;

        this->dt = dt;
    }

    TimeIntegrator::~TimeIntegrator()
    {
    }

    TimeIntegrator *CreateIntegrator(assembly::MonoFieldAssembler *assembler, Float dt, int type)
    {
        TimeIntegrator *integrator = nullptr;

        switch (type)
        {
        case IMPLICIT_EULER:
        {
            integrator = new ImplicitEuler(assembler, dt);
            break;
        }

        default:
            break;
        }

        return integrator;
    }

}