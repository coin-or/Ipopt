// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

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
#include "IpZeroSymMatrix.hpp"
#include "IpCompoundMatrix.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpTransposeMatrix.hpp"
#include "IpExpandedMultiVectorMatrix.hpp"

#include "IpDenseVector.hpp"
#include "IpCompoundVector.hpp"

#include "IpBlas.hpp"

namespace Ipopt
{

  Index TripletHelper::GetNumberEntries(const Matrix& matrix)
  {
    const Matrix* mptr = &matrix;
    const GenTMatrix* gent = dynamic_cast<const GenTMatrix*>(mptr);
    if (gent) {
      return gent->Nonzeros();
    }

    const SymTMatrix* symt = dynamic_cast<const SymTMatrix*>(mptr);
    if (symt) {
      return symt->Nonzeros();
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
      return diag->Dim();
    }

    const IdentityMatrix* ident = dynamic_cast<const IdentityMatrix*>(mptr);
    if (ident) {
      return ident->Dim();
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

    const ZeroSymMatrix* zerosym = dynamic_cast<const ZeroSymMatrix*>(mptr);
    if (zerosym) {
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
      return GetNumberEntries_(*trans);
    }

    const ExpandedMultiVectorMatrix* expmv = dynamic_cast<const ExpandedMultiVectorMatrix*>(mptr);
    if (expmv) {
      return GetNumberEntries_(*expmv);
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to TripletHelper::GetNumberEntries");
  }

  void TripletHelper::FillRowCol(Index n_entries, const Matrix& matrix, Index* iRow, Index* jCol, Index row_offset/*=0*/, Index col_offset/*=0*/)
  {
    const Matrix* mptr = &matrix;
    const GenTMatrix* gent = dynamic_cast<const GenTMatrix*>(mptr);
    if (gent) {
      FillRowCol_(n_entries, *gent, row_offset, col_offset, iRow, jCol);
      return;
    }

    const SymTMatrix* symt = dynamic_cast<const SymTMatrix*>(mptr);
    if (symt) {
      FillRowCol_(n_entries, *symt, row_offset, col_offset, iRow, jCol);
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

    const ZeroSymMatrix* zerosym = dynamic_cast<const ZeroSymMatrix*>(mptr);
    if (zerosym) {
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

    const ExpandedMultiVectorMatrix* expmv = dynamic_cast<const ExpandedMultiVectorMatrix*>(mptr);
    if (expmv) {
      FillRowCol_(n_entries, *expmv, row_offset, col_offset, iRow, jCol);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to TripletHelper::FillRowCol");
  }

  void TripletHelper::FillValues(Index n_entries, const Matrix& matrix, Number* values)
  {
    const Matrix* mptr = &matrix;
    const GenTMatrix* gent = dynamic_cast<const GenTMatrix*>(mptr);
    if (gent) {
      FillValues_(n_entries, *gent, values);
      return;
    }

    const SymTMatrix* symt = dynamic_cast<const SymTMatrix*>(mptr);
    if (symt) {
      FillValues_(n_entries, *symt, values);
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

    const ZeroSymMatrix* zerosym = dynamic_cast<const ZeroSymMatrix*>(mptr);
    if (zerosym) {
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

    const ExpandedMultiVectorMatrix* expmv = dynamic_cast<const ExpandedMultiVectorMatrix*>(mptr);
    if (expmv) {
      FillValues_(n_entries, *expmv, values);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_MATRIX_TYPE,"Unknown matrix type passed to TripletHelper::FillValues");
  }

  Index TripletHelper::GetNumberEntries_(const SumMatrix& matrix)
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

  Index TripletHelper::GetNumberEntries_(const SumSymMatrix& matrix)
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

  Index TripletHelper::GetNumberEntries_(const CompoundMatrix& matrix)
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

  Index TripletHelper::GetNumberEntries_(const CompoundSymMatrix& matrix)
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


  Index TripletHelper::GetNumberEntries_(const TransposeMatrix& matrix)
  {
    return GetNumberEntries(*matrix.OrigMatrix());
  }

  Index TripletHelper::GetNumberEntries_(const ExpandedMultiVectorMatrix& matrix)
  {
    Index nRows = matrix.NRows();
    Index dimVec = matrix.ExpandedMultiVectorMatrixOwnerSpace()->RowVectorSpace()->Dim();
    return nRows*dimVec;
  }

  void TripletHelper::FillRowCol_(Index n_entries, const GenTMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.Nonzeros());
    const Index* irow = matrix.Irows();
    const Index* jcol = matrix.Jcols();
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = irow[i] + row_offset;
      jCol[i] = jcol[i] + col_offset;
    }
  }

  void TripletHelper::FillValues_(Index n_entries, const GenTMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.Nonzeros());
    const Number* vals = matrix.Values();
    for (Index i=0; i<n_entries; i++) {
      values[i] = vals[i];
    }
  }

  void TripletHelper::FillRowCol_(Index n_entries, const SymTMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.Nonzeros());
    const Index* irow = matrix.Irows();
    const Index* jcol = matrix.Jcols();
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = irow[i] + row_offset;
      jCol[i] = jcol[i] + col_offset;
    }
  }

  void TripletHelper::FillValues_(Index n_entries, const SymTMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.Nonzeros());
    matrix.FillValues(values);
  }

  void TripletHelper::FillRowCol_(Index n_entries, const DiagMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.Dim());
    row_offset++;
    col_offset++;
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = i + row_offset;
      jCol[i] = i + col_offset;
    }
  }

  void TripletHelper::FillValues_(Index n_entries, const DiagMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.Dim());
    SmartPtr<const Vector> v = matrix.GetDiag();
    FillValuesFromVector(n_entries, *v, values);
  }

  void TripletHelper::FillRowCol_(Index n_entries, const IdentityMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    DBG_ASSERT(n_entries == matrix.Dim());
    row_offset++;
    col_offset++;
    for (Index i=0; i<n_entries; i++) {
      iRow[i] = i + row_offset;
      jCol[i] = i + col_offset;
    }
  }

  void TripletHelper::FillValues_(Index n_entries, const IdentityMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.Dim());
    Number factor = matrix.GetFactor();
    for (Index i=0; i<n_entries; i++) {
      values[i] = factor;
    }
  }

  void TripletHelper::FillRowCol_(Index n_entries, const ExpansionMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
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

  void TripletHelper::FillValues_(Index n_entries, const ExpansionMatrix& matrix, Number* values)
  {
    DBG_ASSERT(n_entries == matrix.NCols());
    for (Index i=0; i<n_entries; i++) {
      values[i] = 1.0;
    }
  }

  void TripletHelper::FillRowCol_(Index n_entries, const SumMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
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

  void TripletHelper::FillValues_(Index n_entries, const SumMatrix& matrix, Number* values)
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

  void TripletHelper::FillRowCol_(Index n_entries, const SumSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
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

  void TripletHelper::FillValues_(Index n_entries, const SumSymMatrix& matrix, Number* values)
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

  void TripletHelper::FillRowCol_(Index n_entries, const CompoundMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
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

  void TripletHelper::FillValues_(Index n_entries, const CompoundMatrix& matrix, Number* values)
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

  void TripletHelper::FillRowCol_(Index n_entries, const CompoundSymMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
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

  void TripletHelper::FillValues_(Index n_entries, const CompoundSymMatrix& matrix, Number* values)
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

  void TripletHelper::FillValuesFromVector(Index dim, const Vector& vector, Number* values)
  {
    DBG_ASSERT(dim == vector.Dim());
    const DenseVector* dv = dynamic_cast<const DenseVector*>(&vector);
    if (dv) {
      if (dv->IsHomogeneous()) {
        Number scalar = dv->Scalar();
        IpBlasDcopy(dim, &scalar, 0, values, 1);
      }
      else {
        const Number* dv_vals = dv->Values();
        IpBlasDcopy(dim, dv_vals, 1, values, 1);
      }
      return;
    }

    const CompoundVector* cv = dynamic_cast<const CompoundVector*>(&vector);
    if (cv) {
      Index ncomps = cv->NComps();
      Index total_dim = 0;
      for (Index i=0; i<ncomps; i++) {
        SmartPtr<const Vector> comp = cv->GetComp(i);
        Index comp_dim = comp->Dim();
        FillValuesFromVector(comp_dim, *comp, values);
        values += comp_dim;
        total_dim += comp_dim;
      }
      DBG_ASSERT(total_dim == dim);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"Unknown vector type passed to TripletHelper::FillValues");
  }

  void TripletHelper::FillRowCol_(Index n_entries, const ScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, row_offset, col_offset);
  }

  void TripletHelper::FillValues_(Index n_entries, const ScaledMatrix& matrix, Number* values)
  {
    // ToDo:
    // This method can be made much more efficient for ScaledMatrix with GenTMatrix
    // contained

    // Get the matrix values
    FillValues(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), values);

    // Scale the values
    // To Do : This assumes 1-base values (like the TMatrices)
    Index* iRow = new Index[n_entries];
    Index* jCol = new Index[n_entries];
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, 0, 0);

    if (IsValid(matrix.RowScaling())) {
      Index n_rows = matrix.NRows();
      Number* row_scaling = new Number[n_rows];
      FillValuesFromVector(n_rows, *matrix.RowScaling(), row_scaling);
      for (Index i=0; i<n_entries; i++) {
        values[i] *= row_scaling[iRow[i]-1];
      }
      delete [] row_scaling;
    }

    if (IsValid(matrix.ColumnScaling())) {
      Index n_cols = matrix.NCols();
      Number* col_scaling = new Number[n_cols];
      FillValuesFromVector(n_cols, *matrix.ColumnScaling(), col_scaling);
      for (Index i=0; i<n_entries; i++) {
        values[i] *= col_scaling[jCol[i]-1];
      }
      delete [] col_scaling;
    }

    delete [] iRow;
    delete [] jCol;
  }

  void TripletHelper::FillRowCol_(Index n_entries, const SymScaledMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, row_offset, col_offset);
  }

  void TripletHelper::FillValues_(Index n_entries, const SymScaledMatrix& matrix, Number* values)
  {
    // ToDo:
    // This method can be made much more efficient for ScaledMatrix with SymTMatrix
    // contained

    // Get the matrix values
    FillValues(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), values);

    // Scale the values
    // To Do : This assumes 1-base values (like the TMatrices)
    Index* iRow = new Index[n_entries];
    Index* jCol = new Index[n_entries];
    FillRowCol(n_entries, *GetRawPtr(matrix.GetUnscaledMatrix()), iRow, jCol, 0, 0);

    if (IsValid(matrix.RowColScaling())) {
      Index n_dim = matrix.NRows();
      Number* scaling = new Number[n_dim];
      FillValuesFromVector(n_dim, *matrix.RowColScaling(), scaling);
      for (Index i=0; i<n_entries; i++) {
        values[i] *= scaling[iRow[i]-1];
        values[i] *= scaling[jCol[i]-1];
      }
      delete [] scaling;
    }

    delete [] iRow;
    delete [] jCol;
  }

  void TripletHelper::FillRowCol_(Index n_entries, const TransposeMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    FillRowCol(n_entries, *matrix.OrigMatrix(), jCol, iRow,
               col_offset, row_offset);
  }

  void TripletHelper::FillValues_(Index n_entries, const TransposeMatrix& matrix, Number* values)
  {
    FillValues(n_entries, *matrix.OrigMatrix(), values);
  }

  void TripletHelper::FillRowCol_(Index n_entries, const ExpandedMultiVectorMatrix& matrix, Index row_offset, Index col_offset, Index* iRow, Index* jCol)
  {
    row_offset++;
    col_offset++;
    const Index nRows = matrix.NRows();
    SmartPtr<const ExpansionMatrix> P = matrix.GetExpansionMatrix();
    if (IsValid(P)) {
      const Index* exppos = P->ExpandedPosIndices();
      const Index nExp = P->NCols();
      DBG_ASSERT(n_entries = nRows*nExp);
      for (Index irow = row_offset; irow < row_offset+nRows; irow++) {
        for (Index j=0; j<nExp; j++) {
          *(iRow++) = irow;
          *(jCol++) = col_offset + exppos[j];
        }
      }
    }
    else {
      const Index nCols =  matrix.NCols();
      DBG_ASSERT(n_entries = nRows*nCols);
      for (Index irow = row_offset; irow < row_offset+nRows; irow++) {
        for (Index jcol = col_offset; jcol < col_offset+nCols; jcol++) {
          *(iRow++) = irow;
          *(jCol++) = jcol;
        }
      }
    }
  }

  void TripletHelper::FillValues_(Index n_entries, const ExpandedMultiVectorMatrix& matrix, Number* values)
  {
    const Index nRows = matrix.NRows();
    SmartPtr<const ExpansionMatrix> P = matrix.GetExpansionMatrix();
    if (IsValid(P)) {
      const Index nExp = P->NCols();
      double* vecvals = new double[nExp];
      for (Index i=0; i<nRows; i++) {
        SmartPtr<const Vector> vec = matrix.GetVector(i);
        if (IsValid(vec)) {
          DBG_ASSERT(vec->Dim() == nExp);
          FillValuesFromVector(nExp, *vec, vecvals);
          for (Index j=0; j<nExp; j++) {
            *(values++) = vecvals[j];
          }
        }
        else {
          for (Index j=0; j<nExp; j++) {
            *(values++) = 0.;
          }
        }
      }
      delete [] vecvals;
    }
    else {
      const Index nCols =  matrix.NCols();
      for (Index i=0; i<nRows; i++) {
        SmartPtr<const Vector> vec = matrix.GetVector(i);
        if (IsValid(vec)) {
          DBG_ASSERT(vec->Dim() == nCols);
          FillValuesFromVector(nCols, *vec, values);
          values += nCols;
        }
        else {
          for (Index j=0; j<nCols; j++) {
            *(values++) = 0.;
          }
        }
      }
    }
  }

  void TripletHelper::PutValuesInVector(Index dim, const double* values, Vector& vector)
  {
    DBG_ASSERT(dim == vector.Dim());
    DenseVector* dv = dynamic_cast<DenseVector*>(&vector);
    if (dv) {
      Number* dv_vals = dv->Values();
      IpBlasDcopy(dim, values, 1, dv_vals, 1);
      return;
    }

    CompoundVector* cv = dynamic_cast<CompoundVector*>(&vector);
    if (cv) {
      Index ncomps = cv->NComps();
      Index total_dim = 0;
      for (Index i=0; i<ncomps; i++) {
        SmartPtr<Vector> comp = cv->GetCompNonConst(i);
        Index comp_dim = comp->Dim();
        PutValuesInVector(comp_dim, values, *comp);
        values += comp_dim;
        total_dim += comp_dim;
      }
      DBG_ASSERT(total_dim == dim);
      return;
    }

    THROW_EXCEPTION(UNKNOWN_VECTOR_TYPE,"Unknown vector type passed to TripletHelper::PutValuesInVector");
  }

} // namespace Ipopt

