#include "petsc.h"

namespace sfem::la
{
  //=============================================================================
  Vector::Vector(int n_local, int n_global, const std::vector<int> &ghosts)
  {
    VecCreateGhost(SFEM_COMM_WORLD, n_local, n_global, ghosts.size(),
                   ghosts.data(), &_x);
    SetAll(0.0);
  }
  //=============================================================================
  Vector::Vector(Vec x, bool inc_ref_count) : _x(x)
  {
    if (!_x)
    {
      // ERROR
    }

    if (inc_ref_count)
    {
      PetscObjectReference((PetscObject)_x);
    }
  }
  //=============================================================================
  Vector::Vector(Vector &&v)
  {
    _x = v._x;
    v._x = nullptr;
  }
  //=============================================================================
  Vector &Vector::operator=(Vector &&v)
  {
    if (this != &v)
    {
      if (_x)
      {
        VecDestroy(&_x);
      }
      _x = v._x;
      v._x = nullptr;
    }

    return *this;
  }
  //=============================================================================
  Vector::~Vector()
  {
    if (_x)
    {
      VecDestroy(&_x);
    }
  }
  //=============================================================================
  int Vector::GetLocalSize() const
  {
    PetscInt n;
    VecGetLocalSize(_x, &n);
    return n;
  }
  //=============================================================================
  int Vector::GetGlobalSize() const
  {
    PetscInt n;
    VecGetSize(_x, &n);
    return n;
  }
  //=============================================================================
  Vec Vector::Get() const { return _x; }
  //=============================================================================
  Vector Vector::Copy() const
  {
    Vec _y;
    VecDuplicate(_x, &_y);
    VecCopy(_x, _y);
    Vector v(_y, false);
    return v;
  }
  //=============================================================================
  void Vector::SetAll(Float value) { VecSet(_x, value); }
  //=============================================================================
  void Vector::AddValues(int n, const int idxs[], const Float values[])
  {
    VecSetValues(_x, n, idxs, values, ADD_VALUES);
  }
  //=============================================================================
  void Vector::InsertValues(int n, const int idxs[], const Float values[])
  {
    VecSetValues(_x, n, idxs, values, INSERT_VALUES);
  }
  //=============================================================================
  void Vector::Assemble()
  {
    VecAssemblyBegin(_x);
    VecAssemblyEnd(_x);
  }
  //=============================================================================
  std::vector<Float> Vector::GetLocalValues() const
  {
    Vec x_local;
    PetscInt n;
    Float *_values;
    std::vector<Float> values;

    VecGhostUpdateBegin(_x, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(_x, INSERT_VALUES, SCATTER_FORWARD);

    VecGhostGetLocalForm(_x, &x_local);
    VecGetLocalSize(x_local, &n);
    VecGetArray(x_local, &_values);

    values.resize(n);
    for (auto i = 0; i < n; i++)
    {
      values[i] = _values[i];
    }

    VecRestoreArray(x_local, &_values);
    VecGhostRestoreLocalForm(_x, &x_local);

    return values;
  }
  //=============================================================================
  Matrix::Matrix(int n_local, int n_global, const std::vector<int> &diag_nnz,
                 const std::vector<int> &off_diag_nnz)
  {
    MatCreateAIJ(SFEM_COMM_WORLD, n_local, n_local, n_global, n_global,
                 PETSC_DECIDE, diag_nnz.data(), PETSC_DECIDE, nullptr,
                 &_matA); // FIXME: off_diag_nnz error
    MatSetFromOptions(_matA);
  }
  //=============================================================================
  Matrix::Matrix(Mat A, bool inc_ref_count)
      : _matA(A)
  {
    if (!_matA)
    {
      // ERROR
    }

    if (inc_ref_count)
    {
      PetscObjectReference((PetscObject)_matA);
    }
  }

  //=============================================================================
  Matrix::Matrix(Matrix &&A)
  {
    _matA = A._matA;
    A._matA = nullptr;
  }
  //=============================================================================
  Matrix &Matrix::operator=(Matrix &&A)
  {
    if (this != &A)
    {
      if (_matA)
      {
        MatDestroy(&_matA);
      }
      _matA = A._matA;
      A._matA = nullptr;
    }

    return *this;
  }
  //=============================================================================
  Matrix::~Matrix()
  {
    if (_matA)
    {
      MatDestroy(&_matA);
    }
  }
  //=============================================================================
  int Matrix::GetLocalSize() const
  {
    PetscInt n;
    MatGetLocalSize(_matA, &n, nullptr);
    return n;
  }
  //=============================================================================
  int Matrix::GetGlobalSize() const
  {
    PetscInt n;
    MatGetSize(_matA, &n, nullptr);
    return n;
  }
  //=============================================================================
  Mat Matrix::Get() const { return _matA; }
  //=============================================================================
  void Matrix::AddValues(int n, const int idxs[], const Float values[])
  {
    MatSetValues(_matA, n, idxs, n, idxs, values, ADD_VALUES);
  }
  //=============================================================================
  void Matrix::Assemble()
  {
    MatAssemblyBegin(_matA, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(_matA, MAT_FINAL_ASSEMBLY);
  }
  //=============================================================================
  void Matrix::ZeroRowsColumns(int n, const int idxs[], const Float values[],
                               Vector *x, Vector *b)
  {
    x->InsertValues(n, idxs, values);
    MatZeroRowsColumns(_matA, n, idxs, 1.0, x->Get(), b->Get());
  }
  //=============================================================================
  void Matrix::Reset() { MatResetPreallocation(_matA); }
  //=============================================================================
  void PetscSolve(Matrix *A, Vector *b, Vector *x)
  {
    KSP solver;
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, A->Get(), A->Get());
    KSPSetFromOptions(solver);
    KSPSolve(solver, b->Get(), x->Get());
    KSPDestroy(&solver);
  }

} // namespace sfem::la