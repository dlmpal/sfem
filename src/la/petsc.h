#pragma once

#include "../common/config.h"
#include "../third_party/sfem_petsc.h"

namespace sfem::la
{
    /// @brief Thin wrapper around the PETSc Vec
    class Vector
    {
    public:
        /// @brief Create a Vector
        /// @param n_local Local size
        /// @param n_global Global size
        /// @param ghosts Ghost indeces
        Vector(int n_local, int n_global, const std::vector<int> &ghosts);

        /// @brief Create a Vector from an existing PETSc Vec
        /// @param x Existing PETSc Vec
        /// @param inc_ref_count Whether to increase the reference count for x
        Vector(Vec x, bool inc_ref_count);

        // Prevent unnecessary copies
        Vector(const Vector &) = delete;
        Vector &operator=(const Vector &) = delete;

        /// @brief Move constructor
        Vector(Vector &&);

        /// @brief Move assignment
        Vector &operator=(Vector &&);

        /// @brief Destructor
        ~Vector();

        /// @brief Get the local size
        int GetLocalSize() const;

        /// @brief Get the global size
        int GetGlobalSize() const;

        /// @brief Get the underlying Vec
        Vec Get() const;

        /// @brief Copy the vector
        Vector Copy() const;

        /// @brief Set all vector values
        void SetAll(Float value);

        /// @brief Add values into the vector
        /// @param n Size
        /// @param idxs Indeces
        /// @param values Values
        void AddValues(int n, const int idxs[], const Float values[]);

        /// @brief Insert values into the vector, overriding existing ones
        /// @param n Size
        /// @param idxs Indeces
        /// @param values Values
        void InsertValues(int n, const int idxs[], const Float values[]);

        /// @brief Assemble the vector
        void Assemble();

        /// @brief Get the local values
        /// @note The values for ghost indeces are included
        std::vector<Float> GetLocalValues() const;

    private:
        /// @brief Underlying PETSc Vec
        Vec _x;
    };

    /// @brief Thin wrapper around the PETSc Mat
    class Matrix
    {
    public:
        /// @brief Create a Matrix
        /// @param n_local Local size
        /// @param n_global Global size
        /// @param diag_nnz Number of diagonal non-zeros
        /// @param off_diag_nnz Number of off-diagonal non-zeros
        Matrix(int n_local, int n_global, const std::vector<int> &diag_nnz, const std::vector<int> &off_diag_nnz);

        /// @brief Create a Matrix from an existing PETSc Mat
        /// @param A Existing PETSc Mat
        /// @param inc_ref_count Whether to increase the ref count for A
        Matrix(Mat A, bool inc_ref_count);

        // Avoid unnecessary copies
        Matrix(const Matrix &) = delete;
        Matrix &operator=(const Matrix &) = delete;

        /// @brief Move constructor
        Matrix(Matrix &&);

        /// @brief Move assignment
        Matrix &operator=(Matrix &&);

        /// @brief Destructor
        ~Matrix();

        /// @brief Get the local size
        int GetLocalSize() const;

        /// @brief Get the global size
        int GetGlobalSize() const;

        /// @brief Get the underlying PETSc Mat
        Mat Get() const;

        /// @brief Reset the preallocation
        /// @note Call before re-assembling
        void Reset();

        /// @brief Add values into the Matrix
        /// @param n Size
        /// @param idxs Indeces
        /// @param values Values
        void AddValues(int n, const int idxs[], const Float values[]);

        /// @brief Assemble the Matrix
        void Assemble();

        /// @brief Remove rows and columns corresponding to fixed DOF
        /// @param n Size
        /// @param idxs Indeces
        /// @param values Values
        /// @param x Solution vector
        /// @param b RHS vector
        void ZeroRowsColumns(int n, const int idxs[], const Float values[], Vector *x, Vector *b);

    private:
        /// @brief Underlying PETSc Mat
        Mat _matA;
    };

    /// @brief Solve the linear system Ax = b using PETSc's KSP solvers
    void PetscSolve(Matrix *A, Vector *b, Vector *x);
}