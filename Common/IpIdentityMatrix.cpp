// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIdentityMatrix.hpp"

namespace Ipopt
{

  IdentityMatrix::IdentityMatrix(const SymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      factor_(1.0)
  {}

  IdentityMatrix::~IdentityMatrix()
  {}

  Index IdentityMatrix::Dim() const
  {
    DBG_ASSERT(NRows() == NCols());
    return NRows();
  }

  void IdentityMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NRows() == NCols());
    DBG_ASSERT(NRows()==x.Dim());
    DBG_ASSERT(NCols()==y.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    y.Axpy(alpha*factor_, x);
  }

  void IdentityMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    DBG_ASSERT(NRows() == NCols());
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sIdentityMatrix \"%s\" with %d rows and columns and the factor %23.16e.\n",
            prefix.c_str(), name.c_str(), NRows(), factor_);
  }
} // namespace Ipopt
