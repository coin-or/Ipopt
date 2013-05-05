// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter           IBM    2008-08-25

#include "IpTransposeMatrix.hpp"

namespace Ipopt
{

  TransposeMatrix::TransposeMatrix(const TransposeMatrixSpace* owner_space)
      :
      Matrix(owner_space)
  {
    orig_matrix_ = owner_space->MakeNewOrigMatrix();
  }

  void TransposeMatrix::PrintImpl(const Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category,
                                  const std::string& name,
                                  Index indent,
                                  const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sTransposeMatrix \"%s\" of the following matrix\n",
                         prefix.c_str(), name.c_str());
    std::string new_name = name+"^T";
    orig_matrix_->Print(&jnlst, level, category, new_name,
                        indent+1, prefix);
  }

} // namespace Ipopt
