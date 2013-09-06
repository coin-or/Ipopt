// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpZeroSymMatrix.cpp 2269 2013-05-05 11:32:40Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpZeroSymMatrix.hpp"

namespace Ipopt
{

  ZeroSymMatrix::ZeroSymMatrix(const SymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space)
  {}

  ZeroSymMatrix::~ZeroSymMatrix()
  {}

  void ZeroSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                  Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());

    // Take care of the y part of the addition
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }
  }

  void ZeroSymMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
                                       Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==y.Dim());
    DBG_ASSERT(Dim()==x.Dim());

    // Take care of the y part of the addition
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }
  }

  void ZeroSymMatrix::PrintImpl(const Journalist& jnlst,
                             EJournalLevel level,
                             EJournalCategory category,
                             const std::string& name,
                             Index indent,
                             const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sZeroSymMatrix \"%s\" with %d row and %d column components:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols());
  }
} // namespace Ipopt
