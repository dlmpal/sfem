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

        /// @brief Base constructor
        Kernel(KernelType type);

        /// @brief Destructor
        ~Kernel();

        /// @brief
        KernelType GetType() const;

        /// @brief Get the size of the Kernel output
        int GetSize() const;

        /// @brief Integrate the Kernel
        /// @param[in] fe The FiniteElement with which to integrate
        /// @param[out] kloc Where to place kernel output
        virtual void Integrate(fe::FiniteElement *fe, Float kloc[]);

    protected:
        /// @brief Evaluate the Kernel at the current integration point
        virtual void operator()(Float kloc[]) = 0;

        /// @brief Type
        KernelType type;

        /// @brief FiniteElement with which the Kernel is integrated
        fe::FiniteElement *fe;
    };
}