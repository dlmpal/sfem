// Convert a mesh in native format to VTK.
// The program is executed, inside the directory where the field values are located, as follows:
//   sfemToVTK $mesh_path $time_steps $field_name_1 $field_name_2 ...
// It is expected that for each field and time step there exists a separate file,
// e.g. ($field_name_1)_0 for the first field and first time step.
// The program produces as many VTK files as there are timesteps
// The VTK files are placed in the current directory, and named as sfem_0.vtk, sfem_1.vtk, etc

#include "sfem.h"

int main(int argc, char **argv)
{
    sfem::Initialize("SFEM_TO_VTK", nullptr, nullptr);
    std::string mesh_path = argv[1];
    auto mesh = sfem::io::ReadMesh(mesh_path);
    int time_steps = std::atoi(argv[2]);
    int n_fields = argc - 3;
    std::vector<sfem::field::Field> fields;

    // Initialize all fields
    for (auto i = 0; i < n_fields; i++)
    {
        std::string name = argv[i + 3];
        std::string field_path = name + "_" + std::to_string(0);
        std::ifstream file(field_path);
        if (!file.is_open())
            sfem::error::InvalidFileNameError(field_path, __FILE__, __LINE__);
        int n_vars;
        file >> n_vars;
        fields.push_back(sfem::field::Field(name, n_vars, &mesh));
    }

    // Create a .vtk file for each timestep
    for (auto time = 0; time < time_steps; time++)
    {
        // Update field values every timestep
        for (auto &field : fields)
        {
            std::string field_path = field.GetName() + "_" + std::to_string(time);
            sfem::io::ReadFieldValues(field_path, field);
        }
        std::string vtk_filepath = "sfem_" + std::to_string(time) + ".vtk";
        sfem::io::WriteVTK(vtk_filepath, mesh, fields);
    }
    sfem::Finalize();
    return 0;
}
