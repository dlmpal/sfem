#include "basis.h"
#include "../../common/error.h"

namespace sfem::fe::basis
{
    //=============================================================================
    Basis::~Basis()
    {
    }
    //=============================================================================
    int Basis::GetDim() const
    {
        return dim;
    }
    //=============================================================================
    int Basis::GetNumNodes() const
    {
        return n_nodes;
    }
    //=============================================================================
    int Basis::GetNumQuadPts() const
    {
        return n_qpts;
    }
    //=============================================================================
    Float Basis::GetQuadWeight(int npt) const
    {
        return qwts[npt];
    }
    //=============================================================================
    const std::array<Float, 3> &Basis::GetQuadPt(int npt) const
    {
        return qpts[npt];
    }
    //=============================================================================
    Basis *CreateBasis(const mesh::Cell &cell)
    {
        Basis *basis = nullptr;
        switch (cell.type)
        {
        case mesh::CellType::POINT:
            basis = new basis::PointBasis();
            break;

        case mesh::CellType::LINE:
            switch (cell.order)
            {
            case 1:
                basis = new basis::L1Basis();
                break;
            case 2:
                basis = new basis::L2Basis();
                break;
            }
            break;

        case mesh::CellType::TRIANGLE:
            switch (cell.order)
            {
            case 1:
                basis = new basis::Tri1Basis();
                break;
            case 2:
                basis = new basis::Tri2Basis();
                break;
            }
            break;

        case mesh::CellType::TET:
            switch (cell.order)
            {
            case 1:
                basis = new basis::Tet1Basis();
                break;
            case 2:
                basis = new basis::Tet2Basis();
                break;
            }
            break;
        }
        if (basis == nullptr)
        {
            error::InvalidCellError(cell.id, (int)cell.type, cell.order, __FILE__, __LINE__);
        }
        return basis;
    }
}