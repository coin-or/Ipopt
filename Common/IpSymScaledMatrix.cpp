// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSymScaledMatrix.hpp"

namespace Ipopt
{

  SymScaledMatrix::SymScaledMatrix(const SymScaledMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space)
  {}


  SymScaledMatrix::~SymScaledMatrix()
  {}

  void SymScaledMatrix::MultVectorImpl(Number alpha, const Vector &x,
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

    if (IsValid(owner_space_->RowColScaling())) {
      tmp_x->ElementWiseMultiply(*owner_space_->RowColScaling());
    }

    matrix_->MultVector(1.0, *tmp_x, 0.0, *tmp_y);

    if (IsValid(owner_space_->RowColScaling())) {
      tmp_y->ElementWiseMultiply(*owner_space_->RowColScaling());
    }

    y.Axpy(1.0, *tmp_y);
  }

  void SymScaledMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sSymScaledMatrix \"%s\" of dimension %d x %d:\n",
            prefix.c_str(), name.c_str(), NRows(), NCols());
    owner_space_->RowColScaling()->Print(fp, name+"_row_col_scaling", indent, prefix);
    if (IsValid(matrix_)) {
      matrix_->Print(fp, name+"_unscaled_matrix", indent, prefix);
    }
    else {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "unscaled matrix is NULL\n");
    }
  }

} // namespace Ipopt
