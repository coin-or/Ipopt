// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-06-17
//             (based on IpTripletHelper.cpp rev 1312)

#include "IpParTripletHelper.hpp"
#include "IpTripletHelper.hpp"

#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"
#include "IpDiagMatrix.hpp"
#include "IpIdentityMatrix.hpp"
#include "IpExpansionMatrix.hpp"
#include "IpScaledMatrix.hpp"
#include "IpSymScaledMatrix.hpp"
#include "IpSumMatrix.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpZeroMatrix.hpp"
#include "IpCompoundMatrix.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpTransposeMatrix.hpp"

#include "IpDenseVector.hpp"
#include "IpCompoundVector.hpp"

#include "IpParVector.hpp"
#include "IpParGenMatrix.hpp"
#include "IpParSymMatrix.hpp"
#include "IpParExpansionMatrix.hpp"

#include "IpBlas.hpp"

// FIXME - proper header files
//extern "C" {
#define MPICH_SKIP_MPICXX
#include "mpi.h"
//}

namespace Ipopt
{

  static void IdentityRange(Index Dim, Index& start, Index& end)
  {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    start = (int)(Dim * ((double)(rank)/(double)(size)));
    end = (int)(Dim * ((double)(rank+1)/(double)(size))) - 1;
  }

  Index ParTripletHelper::GetNumberEntries(const Matrix& matrix)
  {
    const Matrix* mptr = &matrix;
    const ParGenMatrix* pgen = dynamic_cast<const ParGenMatrix*>(mptr);
    if (pgen) {
      return pgen->LocalMatrix()->Nonzeros();
    }

    const ParSymMatrix* psym = dynamic_cast<const ParSymMatrix*>(mptr);
    if (psym) {
      return psym->LocalMatrix()->Nonzeros();
    }

    const ScaledMatrix* scaled = dynamic_cast<const ScaledMatrix*>(mptr);
    if (scaled) {
      return GetNumberEntries(*GetRawPtr(scaled->GetUnscaledMatrix()));
    }

    const SymScaledMatrix* symscaled = dynamic_cast<const SymScaledMatrix*>(mptr);
    if (symscaled) {
      return GetNumberEntries(*GetRawPtr(symscaled->GetUnscaledMatrix()));
    }

    const DiagMatrix* diag = dynamic_cast<const DiagMatrix*>(mptr);
    if (diag) {
      return GetNumberEntries_(*diag->GetDiag());
    }

    const IdentityMatrix* ident = dynamic_cast<const IdentityMatrix*>(mptr);
    if (ident) {
      // TODO: Is this how we want to do this?
      Index start, end;
      IdentityRange(ident->Dim(), start, end);
      return end-start;
    }

    const ParExpansionMatrix* pexp = dynamic_cast<const ParExpansionMatrix*>(mptr);
    if (pexp) {
      return pexp->LocalMatrix()->NCols();
    }

    const ExpansionMatrix* exp = dynamic_cast<const ExpansionMatrix*>(mptr);
    if (exp) {
      return exp->NCols();
    }

    const SumMatrix* sum = dynamic_cast<const SumMatrix*>(mptr);
    if (sum) {
      return GetNumberEntries_(*sum);
    }

    const SumSymMatrix* sumsym = dynamic_cast<const SumSymMatrix*>(mptr);
    if (sumsym) {
      return GetNumberEntries_(*sumsym);
    }

    const ZeroMatrix* zero = dynamic_cast<const ZeroMatrix*>(mptr);
    if (zero) {
      return 0;
    }

    const CompoundMatrix* cmpd = dynamic_cast<const CompoundMatrix*>(mptr);
    if (cmpd) {
      return GetNumberEntries_(*cmpd);
    }

    const CompoundSymMatrix* cmpd_sym = dynamic_cast<const CompoundSymMatrix*>(mptr);
    if (cmpd_sym) {
      return GetNumberEntries_(*cmpd_sym);
    }

    const TransposeMatrix* trans = dynamic_cast<const TransposeMatrix*>(mptr);
    if (trans) {
      return GetNumberEntries(*trans->OrigMatrix());
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to ParTripletHelper::GetNumberEntries");
  }

  void ParTripletHelper::FillRowCol(Index n_entries, const Matrix& matrix, Index* iRow, Index* jCol, Index row_offset/*=0*/, Index col_offset/*=0*/)
  {
    const Matrix* mptr = &matrix;
    const ParGenMatrix* pgen = dynamic_cast<const ParGenMatrix*>(mptr);
    if (pgen) {
      FillRowCol_(n_entries, *pgen, row_offset, col_offset, iRow, jCol);
      return;
    }

    const ParSymMatrix* psym = dynamic_cast<const ParSymMatrix*>(mptr);
    if (psym) {
      FillRowCol_(n_entries, *psym, row_offset, col_offset, iRow, jCol);
      return;
    }

    const ScaledMatrix* scaled = dynamic_cast<const ScaledMatrix*>(mptr);
    if (scaled) {
      FillRowCol_(n_entries, *scaled, row_offset, col_offset, iRow, jCol);
      return;
    }

    const SymScaledMatrix* symscaled = dynamic_cast<const SymScaledMatrix*>(mptr);
    if (symscaled) {
      FillRowCol_(n_entries, *symscaled, row_offset, col_offset, iRow, jCol);
      return;
    }

    const DiagMatrix* diag = dynamic_cast<const DiagMatrix*>(mptr);
    if (diag) {
      FillRowCol_(n_entries, *diag, row_offset, col_offset, iRow, jCol);
      return;
    }

    const IdentityMatrix* ident = dynamic_cast<const IdentityMatrix*>(mptr);
    if (ident) {
      FillRowCol_(n_entries, *ident, row_offset, col_offset, iRow, jCol);
      return;
    }

    const ParExpansionMatrix* pexp = dynamic_cast<const ParExpansionMatrix*>(mptr);
    if (pexp) {
      FillRowCol_(n_entries, *pexp, row_offset, col_offset, iRow, jCol);
      return;
    }

    const ExpansionMatrix* exp = dynamic_cast<const ExpansionMatrix*>(mptr);
    if (exp) {
      FillRowCol_(n_entries, *exp, row_offset, col_offset, iRow, jCol);
      return;
    }

    const SumMatrix* sum = dynamic_cast<const SumMatrix*>(mptr);
    if (sum) {
      FillRowCol_(n_entries, *sum, row_offset, col_offset, iRow, jCol);
      return;
    }

    const SumSymMatrix* sumsym = dynamic_cast<const SumSymMatrix*>(mptr);
    if (sumsym) {
      FillRowCol_(n_entries, *sumsym, row_offset, col_offset, iRow, jCol);
      return;
    }

    const ZeroMatrix* zero = dynamic_cast<const ZeroMatrix*>(mptr);
    if (zero) {
      DBG_ASSERT(n_entries == 0);
      return;
    }

    const CompoundMatrix* cmpd = dynamic_cast<const CompoundMatrix*>(mptr);
    if (cmpd) {
      FillRowCol_(n_entries, *cmpd, row_offset, col_offset, iRow, jCol);
      return;
    }

    const CompoundSymMatrix* cmpd_sym = dynamic_cast<const CompoundSymMatrix*>(mptr);
    if (cmpd_sym) {
      FillRowCol_(n_entries, *cmpd_sym, row_offset, col_offset, iRow, jCol);
      return;
    }

    const TransposeMatrix* trans = dynamic_cast<const TransposeMatrix*>(mptr);
    if (trans) {
      FillRowCol_(n_entries, *trans, row_offset, col_offset, iRow, jCol);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to ParTripletHelper::FillRowCol");
  }

  void ParTripletHelper::FillValues(Index n_entries, const Matrix& matrix, Number* values)
  {
    const Matrix* mptr = &matrix;
    const ParGenMatrix* pgen = dynamic_cast<const ParGenMatrix*>(mptr);
    if (pgen) {
      FillValues_(n_entries, *pgen, values);
      return;
    }

    const ParSymMatrix* psym = dynamic_cast<const ParSymMatrix*>(mptr);
    if (psym) {
      FillValues_(n_entries, *psym, values);
      return;
    }

    const ScaledMatrix* scaled = dynamic_cast<const ScaledMatrix*>(mptr);
    if (scaled) {
      FillValues_(n_entries, *scaled, values);
      return;
    }

    const SymScaledMatrix* symscaled = dynamic_cast<const SymScaledMatrix*>(mptr);
    if (symscaled) {
      FillValues_(n_entries, *symscaled, values);
      return;
    }

    const DiagMatrix* diag = dynamic_cast<const DiagMatrix*>(mptr);
    if (diag) {
      FillValues_(n_entries, *diag, values);
      return;
    }

    const IdentityMatrix* ident = dynamic_cast<const IdentityMatrix*>(mptr);
    if (ident) {
      FillValues_(n_entries, *ident, values);
      return;
    }

    const ParExpansionMatrix* pexp = dynamic_cast<const ParExpansionMatrix*>(mptr);
    if (pexp) {
      FillValues_(n_entries, *pexp, values);
      return;
    }

    const ExpansionMatrix* exp = dynamic_cast<const ExpansionMatrix*>(mptr);
    if (exp) {
      FillValues_(n_entries, *exp, values);
      return;
    }

    const SumMatrix* sum = dynamic_cast<const SumMatrix*>(mptr);
    if (sum) {
      FillValues_(n_entries, *sum, values);
      return;
    }

    const SumSymMatrix* sumsym = dynamic_cast<const SumSymMatrix*>(mptr);
    if (sumsym) {
      FillValues_(n_entries, *sumsym, values);
      return;
    }

    const ZeroMatrix* zero = dynamic_cast<const ZeroMatrix*>(mptr);
    if (zero) {
      DBG_ASSERT(n_entries == 0);
      return;
    }

    const CompoundMatrix* cmpd = dynamic_cast<const CompoundMatrix*>(mptr);
    if (cmpd) {
      FillValues_(n_entries, *cmpd, values);
      return;
    }

    const CompoundSymMatrix* cmpd_sym = dynamic_cast<const CompoundSymMatrix*>(mptr);
    if (cmpd_sym) {
      FillValues_(n_entries, *cmpd_sym, values);
      return;
    }

    const TransposeMatrix* trans = dynamic_cast<const TransposeMatrix*>(mptr);
    if (trans) {
      FillValues_(n_entries, *trans, values);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to ParTripletHelper::FillValues");
  }

  Index ParTripletHelper::GetNumberEntries_(const SumMatrix& matrix)
  {
    Index n_entries = 0;
    Index nterms = matrix.NTerms();
    for (Index i=0; i<nterms; i++) {
      Number dummy;
      SmartPtr<const Matrix> i_mat;
      matrix.GetTerm(i, dummy, i_mat);
      n_entries += GetNumberEntries(*i_mat);
    }
    return n_entries;
  }

  Index ParTripletHelper::GetNumberEntries_(const SumSymMatrix& matrix)
  {
    Index n_entries = 0;
    Index nterms = matrix.NTerms();
    for (Index i=0; i<nterms; i++) {
      Number dummy;
      SmartPtr<const SymMatrix> i_mat;
      matrix.GetTerm(i, dummy, i_mat);
      n_entries += GetNumberEntries(*i_mat);
    }
    return n_entries;
  }

  Index ParTripletHelper::GetNumberEntries_(const CompoundMatrix& matrix)
  {
    Index n_entries = 0;
    Index nrows = matrix.NComps_Rows();
    Index ncols = matrix.NComps_Cols();
    for (Index i=0; i<nrows; i++) {
      for (Index j=0; j<ncols; j++) {
        SmartPtr<const Matrix> comp = matrix.GetComp(i,j);
        if (IsValid(comp)) {
          n_entries += GetNumberEntries(*comp);
        }
      }
    }
    return n_entries;
  }

  Index ParTripletHelper::GetNumberEntries_(const CompoundSymMatrix& matrix)
  {
    Index n_entries = 0;
    Index dim = matrix.NComps_Dim();
    for (Index i=0; i<dim; i++) {
      for (Index j=0; j<=i; j++) {
        SmartPtr<const Matrix> comp = matrix.GetComp(i,j);
        if (IsValid(comp)) {
          n_entries += GetNumberEntries(*comp);
        }
      }
    }
    return n_entries;
  }

  Index ParTripletHelper::GetNumberEntries_(const Vector& vector)
  {
    const Vector* vptr = &vector;
    const ParVector* p_vec = dynamic_cast<const ParVector*>(vptr);
    if (p_vec) {
      return p_vec->LocalSize();
    }

    const CompoundVector* cmpd_vec = dynamic_cast<const CompoundVector*>(vptr);
    if (cmpd_vec) {
      Index n_entires = 0;
      Index n_comps = cmpd_vec->NComps();
      for (int i=0; i<n_comps; i++) {
	n_entires += GetNumberEntries_(*cmpd_vec->GetComp(i));
      }
      return n_entires;
    }

    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"Unknown vector type passed to ParTripletHelper::GetNumberEntries_(const Vector)");
  }


  void ParTripletHelper::FillRowCol_(Index n_entries, const ParGenMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    TripletHelper::FillRowCol(n_entries, *matrix.LocalMatrix(),
			      iRow, jCol, row_offset+matrix.RowStartPos(),
			      col_offset);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const ParGenMatrix& matrix, Number* values)
  {
    TripletHelper::FillValues(n_entries, *matrix.LocalMatrix(), values);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const ParSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    TripletHelper::FillRowCol(n_entries, *matrix.LocalMatrix(),
			      iRow, jCol,  row_offset, col_offset);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const ParSymMatrix& matrix, Number* values)
  {
    TripletHelper::FillValues(n_entries, *matrix.LocalMatrix(), values);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const DiagMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol_(n_entries, *matrix.GetDiag(), row_offset, col_offset, iRow, jCol);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const DiagMatrix& matrix, Number* values)
  {
    FillValues_(n_entries, *matrix.GetDiag(), values);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const Vector& vector, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    const Vector* vptr = &vector;
    const ParVector* p_vec = dynamic_cast<const ParVector*>(vptr);
    if (p_vec) {
      const Index size = p_vec->LocalSize();
      const Index roffset = row_offset + p_vec->StartPos();
      const Index coffset = col_offset + p_vec->StartPos();
      DBG_ASSERT(n_entries == size);
      for (Index i=1; i<=size; i++) {
	iRow[i] = i + roffset;
	jCol[i] = i + coffset;
      }
      return;
    }

    const CompoundVector* cmpd_vec = dynamic_cast<const CompoundVector*>(vptr);
    if (cmpd_vec) {
      Index n_total_entries = 0;
      Index n_comps = cmpd_vec->NComps();
      Index roffset = row_offset;
      Index coffset = col_offset;
      for (int i=0; i<n_comps; i++) {
	SmartPtr<const Vector> blk_vec = cmpd_vec->GetComp(i);
	Index blk_n_entries = GetNumberEntries_(*blk_vec);
	FillRowCol_(blk_n_entries, *blk_vec, roffset, coffset, iRow, jCol);
	roffset += blk_vec->Dim();
	coffset += blk_vec->Dim();
	iRow += blk_n_entries;
	jCol += blk_n_entries;
	n_total_entries += blk_n_entries;
      }
      DBG_ASSERT(n_total_entries == n_entries);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"Unknown vector type passed to ParTripletHelper::FillRowCol_(const Vector)");
  }

  void ParTripletHelper::FillValues_(Index n_entries, const Vector& vector, Number* values)
  {
    const Vector* vptr = &vector;
    const ParVector* p_vec = dynamic_cast<const ParVector*>(vptr);
    if (p_vec) {
      DBG_ASSERT(n_entries == p_vec->LocalSize());
      TripletHelper::FillValuesFromVector(n_entries, *p_vec->LocalVector(), values);
      return;
    }

    const CompoundVector* cmpd_vec = dynamic_cast<const CompoundVector*>(vptr);
    if (cmpd_vec) {
      Index n_total_entries = 0;
      Index n_comps = cmpd_vec->NComps();
      for (int i=0; i<n_comps; i++) {
	SmartPtr<const Vector> blk_vec = cmpd_vec->GetComp(i);
	Index blk_n_entries = GetNumberEntries_(*blk_vec);
	FillValues_(blk_n_entries, *blk_vec, values);
	values += blk_n_entries;
	n_total_entries += blk_n_entries;
      }
      DBG_ASSERT(n_total_entries == n_entries);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"Unknown vector type passed to ParTripletHelper::FillRowCol_(const Vector)");
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const IdentityMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    Index start, end;
    IdentityRange(matrix.Dim(), start, end);

    DBG_ASSERT(n_entries == end-start);
    row_offset += start;
    col_offset += start;
    for (Index i=1; i<=n_entries; i++) {
      iRow[i] = i + row_offset;
      jCol[i] = i + col_offset;
    }
  }

  void ParTripletHelper::FillValues_(Index n_entries, const IdentityMatrix& matrix, Number* values)
  {
    Index start, end;
    IdentityRange(matrix.Dim(), start, end);

    DBG_ASSERT(n_entries == end-start);
    Number factor = matrix.GetFactor();
    for (Index i=0; i<n_entries; i++) {
      values[i] = factor;
    }
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const ParExpansionMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.LocalMatrix()->NCols());
    const Index* exp_pos = matrix.LocalMatrix()->ExpandedPosIndices();
    row_offset += matrix.RowStartPos() + 1;
    col_offset += matrix.ColStartPos() + 1;
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = exp_pos[i] + row_offset;
      jCol[i] = i + col_offset;
    }
  }

  void ParTripletHelper::FillValues_(Index n_entries, const ParExpansionMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.LocalMatrix()->NCols());
    for (Index i=0; i<n_entries; i++) {
      values[i] = 1.0;
    }
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const ExpansionMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.NCols());
    const Index* exp_pos = matrix.ExpandedPosIndices();
    row_offset++;
    col_offset++;
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = exp_pos[i] + row_offset;
      jCol[i] = i + col_offset;
    }
  }

  void ParTripletHelper::FillValues_(Index n_entries, const ExpansionMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.NCols());
    for (Index i=0; i<n_entries; i++) {
      values[i] = 1.0;
    }
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const SumMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    Index total_n_entries = 0;
    for (Index i=0; i<matrix.NTerms(); i++) {
      // Fill the indices for the individual term
      Number retFactor = 0.0;
      SmartPtr<const Matrix> retTerm;
      matrix.GetTerm(i, retFactor, retTerm);
      Index term_n_entries = GetNumberEntries(*retTerm);
      total_n_entries += term_n_entries;
      FillRowCol(term_n_entries, *retTerm, iRow, jCol, row_offset, col_offset);

      // now shift the iRow, jCol pointers for the next term
      iRow += term_n_entries;
      jCol += term_n_entries;
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const SumMatrix& matrix, Number* values)
  {
    Index total_n_entries = 0;
    for (Index i=0; i<matrix.NTerms(); i++) {
      // Fill the values for the individual term
      Number retFactor = 0.0;
      SmartPtr<const Matrix> retTerm;
      matrix.GetTerm(i, retFactor, retTerm);
      Index term_n_entries = GetNumberEntries(*retTerm);
      total_n_entries += term_n_entries;
      FillValues(term_n_entries, *retTerm, values);

      // Now adjust the values based on the factor
      IpBlasDscal(term_n_entries, retFactor, values, 1);

      // now shift the values pointer for the next term
      values += term_n_entries;
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const SumSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    Index total_n_entries = 0;
    for (Index i=0; i<matrix.NTerms(); i++) {
      // Fill the indices for the individual term
      Number retFactor = 0.0;
      SmartPtr<const SymMatrix> retTerm;
      matrix.GetTerm(i, retFactor, retTerm);
      Index term_n_entries = GetNumberEntries(*retTerm);
      total_n_entries += term_n_entries;
      FillRowCol(term_n_entries, *retTerm, iRow, jCol, row_offset, col_offset);

      // now shift the iRow, jCol pointers for the next term
      iRow += term_n_entries;
      jCol += term_n_entries;
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const SumSymMatrix& matrix, Number* values)
  {
    Index total_n_entries = 0;
    for (Index i=0; i<matrix.NTerms(); i++) {
      // Fill the values for the individual term
      Number retFactor = 0.0;
      SmartPtr<const SymMatrix> retTerm;
      matrix.GetTerm(i, retFactor, retTerm);
      Index term_n_entries = GetNumberEntries(*retTerm);
      total_n_entries += term_n_entries;
      if (retFactor!=0.0) {
        FillValues(term_n_entries, *retTerm, values);

        if (retFactor!=1.) {
          // Now adjust the values based on the factor
          IpBlasDscal(term_n_entries, retFactor, values, 1);
        }
      }
      else {
        const Number zero = 0.;
        IpBlasDcopy(term_n_entries, &zero, 0, values, 1);
      }

      // now shift the values pointer for the next term
      values += term_n_entries;
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const CompoundMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    Index total_n_entries = 0;

    const CompoundMatrixSpace* owner_space = static_cast<const CompoundMatrixSpace*>(GetRawPtr(matrix.OwnerSpace()));
    DBG_ASSERT(dynamic_cast<const CompoundMatrixSpace*>(GetRawPtr(matrix.OwnerSpace())));

    Index c_row_offset = row_offset;
    for (Index i=0; i<matrix.NComps_Rows(); i++) {
      Index c_col_offset = col_offset;
      for (Index j=0; j<matrix.NComps_Cols(); j++) {
        // Fill the indices for the individual term
        SmartPtr<const Matrix> blk_mat = matrix.GetComp(i, j);
        if (IsValid(blk_mat)) {
          Index blk_n_entries = GetNumberEntries(*blk_mat);
          total_n_entries += blk_n_entries;
          FillRowCol(blk_n_entries, *blk_mat, iRow, jCol, c_row_offset, c_col_offset);

          // now shift the iRow, jCol pointers for the next term
          iRow += blk_n_entries;
          jCol += blk_n_entries;
        }
        c_col_offset += owner_space->GetBlockCols(j);
      }
      c_row_offset += owner_space->GetBlockRows(i);
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const CompoundMatrix& matrix, Number* values)
  {
    Index total_n_entries = 0;

    for (Index i=0; i<matrix.NComps_Rows(); i++) {
      for (Index j=0; j<matrix.NComps_Cols(); j++) {
        // Fill the indices for the individual term
        SmartPtr<const Matrix> blk_mat = matrix.GetComp(i, j);
        if (IsValid(blk_mat)) {
          Index blk_n_entries = GetNumberEntries(*blk_mat);
          total_n_entries += blk_n_entries;
          FillValues(blk_n_entries, *blk_mat, values);

          // now shift the values pointer for the next term
          values += blk_n_entries;
        }
      }
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const CompoundSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    Index total_n_entries = 0;

    const CompoundSymMatrixSpace* owner_space = static_cast<const CompoundSymMatrixSpace*>(GetRawPtr(matrix.OwnerSpace()));
    DBG_ASSERT(dynamic_cast<const CompoundSymMatrixSpace*>(GetRawPtr(matrix.OwnerSpace())));

    Index c_row_offset = row_offset;
    for (Index i=0; i<matrix.NComps_Dim(); i++) {
      Index c_col_offset = col_offset;
      for (Index j=0; j<=i; j++) {
        // Fill the indices for the individual term
        SmartPtr<const Matrix> blk_mat = matrix.GetComp(i, j);
        if (IsValid(blk_mat)) {
          Index blk_n_entries = GetNumberEntries(*blk_mat);
          total_n_entries += blk_n_entries;
          FillRowCol(blk_n_entries, *blk_mat, iRow, jCol, c_row_offset, c_col_offset);

          // now shift the iRow, jCol pointers for the next term
          iRow += blk_n_entries;
          jCol += blk_n_entries;
        }
        c_col_offset += owner_space->GetBlockDim(j);
      }
      c_row_offset += owner_space->GetBlockDim(i);
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const CompoundSymMatrix& matrix, Number* values)
  {
    Index total_n_entries = 0;

    for (Index i=0; i<matrix.NComps_Dim(); i++) {
      for (Index j=0; j<=i; j++) {
        // Fill the indices for the individual term
        SmartPtr<const Matrix> blk_mat = matrix.GetComp(i, j);
        if (IsValid(blk_mat)) {
          Index blk_n_entries = GetNumberEntries(*blk_mat);
          total_n_entries += blk_n_entries;
          FillValues(blk_n_entries, *blk_mat, values);

          // now shift the iRow, jCol pointers for the next term
          values += blk_n_entries;
        }
      }
    }
    DBG_ASSERT(total_n_entries == n_entries);
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const ScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, row_offset, col_offset);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const ScaledMatrix& matrix, Number* values)
  {
    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"ParTripletHelper::FillValues_ has not been implemented for ScaledMatrix");
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const SymScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, row_offset, col_offset);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const SymScaledMatrix& matrix, Number* values)
  {
    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"ParTripletHelper::FillValues_ has not been implemented for SymScaledMatrix");
  }

  void ParTripletHelper::FillRowCol_(Index n_entries, const TransposeMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *matrix.OrigMatrix(), jCol, iRow,
               col_offset, row_offset);
  }

  void ParTripletHelper::FillValues_(Index n_entries, const TransposeMatrix& matrix, Number* values)
  {
    FillValues(n_entries, *matrix.OrigMatrix(), values);
  }

} // namespace Ipopt

