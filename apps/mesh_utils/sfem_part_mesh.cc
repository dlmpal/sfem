#include "sfem.h"

int main(int argc, char **argv)
{
    sfem::Initialize(&argc, &argv, "SFEM_MESH_PART");
    std::string mesh_path = argv[1];
    int n_parts = std::atoi(argv[2]);
    std::string partitioner_type(argv[3]);
    auto mesh = sfem::io::ReadMesh(mesh_path);
    auto partitioner = sfem::mesh::CreateMeshPartitioner(partitioner_type, &mesh, n_parts);
    partitioner->PartitionMesh();
    partitioner->WritePartition(mesh_path);
    return 0;
}