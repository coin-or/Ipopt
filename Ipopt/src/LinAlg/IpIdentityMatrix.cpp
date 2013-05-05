// Copyright (C) 2004, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
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
    DBG_ASSERT(NRows() == x.Dim());
    DBG_ASSERT(NCols() == y.Dim());

    y.AddOneVector(alpha*factor_, x, beta);
  }

  bool IdentityMatrix::HasValidNumbersImpl() const
  {
    return IsFiniteNumber(factor_);
  }

  void IdentityMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    if (init) {
      rows_norms.Set(1.);
    }
    else {
      SmartPtr<Vector> v = rows_norms.MakeNew();
      v->Set(1.);
      rows_norms.ElementWiseMax(*v);
    }
  }

  void IdentityMatrix::PrintImpl(const Journalist& jnlst,
                                 EJournalLevel level,
                                 EJournalCategory category,
                                 const std::string& name,
                                 Index indent,
                                 const std::string& prefix) const
  {
    DBG_ASSERT(NRows() == NCols());
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sIdentityMatrix \"%s\" with %d rows and columns and the factor %23.16e.\n",
                         prefix.c_str(), name.c_str(), NRows(), factor_);
  }

  void IdentityMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
                                     const Vector& Z, Vector& X) const
  {
    DBG_ASSERT(NRows() == NCols());
    DBG_ASSERT(NRows() == S.Dim());
    DBG_ASSERT(NCols() == Z.Dim());
    DBG_ASSERT(NCols() == X.Dim());

    X.AddVectorQuotient(alpha, Z, S, 1.);
  }
} // namespace Ipopt
