#include "field.h"
#include "../common/error.h"

namespace sfem::io
{
    //=============================================================================
    void ReadFieldValues(const std::string &path, field::Field &field)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            error::InvalidFileNameError(path, __FILE__, __LINE__);
        }

        int n_vars;
        file >> n_vars;
        if (n_vars != field.GetNumVars())
        {
            error::InvalidSizeError(field.GetNumVars(), n_vars, __FILE__, __LINE__);
        }

        auto mesh = field.GetMesh();
        const auto &global_to_local = mesh->GetNodeGlobalToLocalMapping(mesh::IndexingType::GLOBAL);
        std::vector<Float> values(field.GetNumDof());
        for (auto i = 0; i < mesh->GetNumNodesGlobal(); i++)
        {
            Float vals[n_vars];
            for (auto j = 0; j < n_vars; j++)
                file >> vals[j];

            // Only add the values for owned + ghost DOF
            if (global_to_local.count(i) == 1)
            {
                for (auto j = 0; j < n_vars; j++)
                {
                    values[global_to_local.at(i) * n_vars + j] = vals[j];
                }
            }
        }
        field.SetValues(values);
    }
    //=============================================================================
    void WriteFieldValues(const std::string &path, const field::Field &field, bool assemble_global)
    {
        std::string _path = path;
        std::vector<Float> values;
        if (assemble_global)
        {
            values = field::AssembleGlobalValues(field);
        }
        else
        {
            values = field.GetValues();
            if (Logger::GetInstance().n_procs > 1)
            {
                _path += std::to_string(Logger::GetInstance().proc_rank);
            }
        }
        if (Logger::GetInstance().proc_rank == SFEM_ROOT || assemble_global == false)
        {
            std::ofstream file(_path);
            if (!file.is_open())
            {
                error::InvalidFileNameError(_path, __FILE__, __LINE__);
            }
            file << field.GetNumVars() << "\n";
            for (auto i = 0; i < values.size(); i++)
            {
                file << values[i] << "\n";
            }
        }
    }
}