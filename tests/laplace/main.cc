// Solve the laplace equation in 1d, 2d or 3d.

#include "sfem.h"

using namespace sfem;

extern void LaplaceSolver(int dim);

int main(int argc, char **argv)
{
    Initialize(&argc, &argv, "LAPLACE_SOLVER", "log");
    int dim = std::atoi(argv[1]);
    LaplaceSolver(dim);
    Finalize();
    return 0;
}

void LaplaceSolver(int dim)
{
    auto mesh = io::ReadMesh("Mesh3D");
    mesh.Info();
    field::Field phi("Phi", 1, &mesh);

    phi.AddFixedDof("Fixed", 0, 10.);
    phi.AddFixedDof("Free", 0, 100.);

    kernel::basic::Diffusion diffusion_kernel(1.0);
    assembly::MonoFieldAssembler assembler(&phi, assembly::AssemblyType::STATIC);
    assembler.AddKernel("Body", &diffusion_kernel);

    solver::StaticSolver solver(&assembler);
    solver.Run();
    io::WriteVTK("fields/sfem_" + std::to_string(::Logger::GetInstance().GetProcRank()) + ".vtk", mesh, {phi});
}