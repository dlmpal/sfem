#pragma once

#include "../fe/finite_element.h"

/// @brief Finite element kernels
namespace sfem::kernel
{
    /// @brief Abstract Kernel class
    class Kernel
    {
    public:
        /// @brief Where to add/insert kernel output
        enum class KernelType : int
        {
            LHS = 0,
            RHS = 1,
            BOTH = 2,
            MASS = 3
        };

        /// @brief Destructor
        ~Kernel();

        /// @brief
        virtual KernelType GetType() const = 0;

        /// @brief Get the size of the Kernel output
        int GetSize() const;

        /// @brief Integrate the Kernel
        /// @param fe The FiniteElement with which to integrate
        /// @returns The integrated kernel output
        virtual std::vector<Float> Integrate(fe::FiniteElement *fe);

    protected:
        /// @brief Evaluate the Kernel at the current integration point
        /// @param[out] kloc Kernel output
        virtual void Evaluate(std::vector<Float> &kloc) = 0;

        /// @brief Type
        KernelType type;

        /// @brief FiniteElement with which the Kernel is integrated
        fe::FiniteElement *fe;
    };
}