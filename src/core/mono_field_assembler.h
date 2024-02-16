#pragma once

#include "../fe/finite_element.h"
#include "../la/sfem_la.h"
#include "field.h"
#include "kernel.h"

/// @brief Discretized operator matrix assembly
namespace sfem::assembly
{
    /// @brief Assembly types
    enum class AssemblyType : int
    {
        STATIC = 0,
        DYNAMIC = 1
    };

    /// @brief Matrix assembly for single field problems
    class MonoFieldAssembler
    {
    public:
        /// @brief Create a MonoFieldAssembler
        MonoFieldAssembler(field::Field *field, AssemblyType type);

        /// @brief Get assembly type
        AssemblyType GetType() const;

        /// @brief Get a pointer to the Field
        field::Field *GetField() const;

        /// @brief Get a pointer to the solution vector
        la::Vector *GetSolutionVector() const;

        /// @brief Get a pointer to the RHS or load vector
        la::Vector *GetLoadVector() const;

        /// @brief Get a pointer to the stiffness matrix
        la::Matrix *GetStiffnessMatrix() const;

        /// @brief Get a pointer to the mass matrix
        /// @note For static problems returns a nullptr
        la::Matrix *GetMassMatrix() const;

        /// @brief Add a kernel to a region
        void AddKernel(std::string region_name, kernel::Kernel *kernel);

        /// @brief Assemble the linear system
        void AssembleSystem();

        /// @brief Update the field values from the solution vector
        void UpdateFieldValues();

    private:
        /// @brief The type of assembly
        AssemblyType type;

        /// @brief The Field for which the PDE is solved
        field::Field *field;

        /// @brief The solution vector
        std::unique_ptr<la::Vector> U;

        /// @brief The right-hand-side (load vector) of the linear system
        std::unique_ptr<la::Vector> F;

        /// @brief The left-hand-side (stiffness matrix) of the linear system
        std::unique_ptr<la::Matrix> K;

        /// @brief The mass matrix
        /// @note Not assembled for static problems
        std::unique_ptr<la::Matrix> M;

        /// @brief Map mesh regions to kernels
        /// @note A single region may have more than one kernels defined
        std::map<std::string, std::vector<kernel::Kernel *>> kernel_map;
    };

}