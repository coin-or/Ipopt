// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpDiagMatrix.hpp"

namespace Ipopt
{

  DiagMatrix::DiagMatrix(const SymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space)
  {}

  DiagMatrix::~DiagMatrix()
  {}

  void DiagMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                  Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());
    DBG_ASSERT(IsValid(diag_));

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    SmartPtr<Vector> tmp_vec = y.MakeNew();
    tmp_vec->Copy(x);
    tmp_vec->ElementWiseMultiply(*diag_);
    y.Axpy(alpha, *tmp_vec);
  }

  void DiagMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sDiagMatrix \"%s\" with %d rows and columns, and with diagonal elements:\n",
            prefix.c_str(), name.c_str(), Dim());
    if (IsValid(diag_)) {
      diag_->Print(fp, name, indent, prefix);
    }
    else {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "%sDiagonal elements not set!\n", prefix.c_str());
    }
  }
} // namespace Ipopt
