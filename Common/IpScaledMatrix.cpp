// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpScaledMatrix.hpp"

namespace Ipopt
{

  ScaledMatrix::ScaledMatrix(const ScaledMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {}


  ScaledMatrix::~ScaledMatrix()
  {}

  void ScaledMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                    Number beta, Vector &y) const
  {
    DBG_ASSERT(IsValid(matrix_));

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // need some temporary vectors
    SmartPtr<Vector> tmp_x = x.MakeNewCopy();
    SmartPtr<Vector> tmp_y = y.MakeNew();

    if (IsValid(owner_space_->ColumnScaling())) {
      tmp_x->ElementWiseMultiply(*owner_space_->ColumnScaling());
    }

    matrix_->MultVector(1.0, *tmp_x, 0.0, *tmp_y);

    if (IsValid(owner_space_->RowScaling())) {
      tmp_y->ElementWiseMultiply(*owner_space_->RowScaling());
    }

    y.Axpy(1.0, *tmp_y);
  }

  void ScaledMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
                                         Number beta, Vector &y) const
  {
    DBG_ASSERT(IsValid(matrix_));

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // need some temporary vectors
    SmartPtr<Vector> tmp_x = x.MakeNewCopy();
    SmartPtr<Vector> tmp_y = y.MakeNew();

    if (IsValid(owner_space_->RowScaling())) {
      tmp_x->ElementWiseMultiply(*owner_space_->RowScaling());
    }

    matrix_->TransMultVector(1.0, *tmp_x, 0.0, *tmp_y);

    if (IsValid(owner_space_->ColumnScaling())) {
      tmp_y->ElementWiseMultiply(*owner_space_->ColumnScaling());
    }

    y.Axpy(1.0, *tmp_y);
  }

  void ScaledMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sScaledMatrix \"%s\" of dimension %d x %d:\n",
            prefix.c_str(), name.c_str(), NRows(), NCols());
    if (IsValid(owner_space_->RowScaling())) {
      owner_space_->RowScaling()->Print(fp, name+"_row_scaling", indent, prefix);
    }
    else {
      fprintf(fp, "RowScaling is NULL\n");
    }
    if (IsValid(matrix_)) {
      matrix_->Print(fp, name+"_unscaled_matrix", indent, prefix);
    }
    else {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "unscaled matrix is NULL\n");
    }
    if (IsValid(owner_space_->ColumnScaling())) {
      owner_space_->ColumnScaling()->Print(fp, name+"_column_scaling", indent, prefix);
    }
    else {
      fprintf(fp, "ColumnScaling is NULL\n");
    }
  }

  void ScaledMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
                                   const Vector& Z, Vector& X) const
  {
    DBG_ASSERT(false && "Got the ScaledMatrix::AddMSinvZImpl.  Should implement specialized method!");

    SmartPtr<Vector> tmp = S.MakeNew();
    tmp->AddVectorQuotient(1., Z, S, 0.);
    MultVector(alpha, *tmp, 1., X);
  }

  void ScaledMatrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
                                        const Vector& R, const Vector& Z,
                                        const Vector& D, Vector& X) const
  {
    DBG_ASSERT(false && "Got the ScaledMatrix::SinvBlrmZMTdBrImpl.  Should implement specialized method!");

    TransMultVector(alpha, D, 0., X);
    X.ElementWiseMultiply(Z);
    X.Axpy(1., R);
    X.ElementWiseDivide(S);
  }


  ScaledMatrixSpace::ScaledMatrixSpace(
    const SmartPtr<const Vector>& row_scaling,
    bool row_scaling_reciprocal,
    const SmartPtr<const MatrixSpace>& unscaled_matrix_space,
    const SmartPtr<const Vector>& column_scaling,
    bool column_scaling_reciprocal)
      :
      MatrixSpace(unscaled_matrix_space->NRows(),
                  unscaled_matrix_space->NCols()),
      unscaled_matrix_space_(unscaled_matrix_space)
  {
    if (IsValid(row_scaling)) {
      row_scaling_ = row_scaling->MakeNewCopy();
      if (row_scaling_reciprocal) {
        row_scaling_->ElementWiseReciprocal();
      }
    }
    else {
      row_scaling_ = NULL;
    }

    if (IsValid(column_scaling)) {
      column_scaling_ = column_scaling->MakeNewCopy();
      if (column_scaling_reciprocal) {
        column_scaling_->ElementWiseReciprocal();
      }
    }
    else {
      column_scaling_ = NULL;
    }
  }
} // namespace Ipopt
