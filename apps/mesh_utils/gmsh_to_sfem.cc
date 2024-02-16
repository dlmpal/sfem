// Convert a gmsh file to native sfem mesh format
// .msh2 files are suported
// The program is executed as follows:
//   gmshToSfem $gmsh_path $mesh_path

#include "sfem.h"

int main(int argc, char **argv)
{
    sfem::Initialize("GMSH_TO_SFEM", nullptr, nullptr);
    std::string gmsh_path = argv[1];
    std::string mesh_path = argv[2];
    auto mesh = sfem::io::ReadGmsh(gmsh_path);
    sfem::io::WriteMesh(mesh_path, mesh);
    sfem::Finalize();
    return 0;
}
