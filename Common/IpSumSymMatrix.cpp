// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSumSymMatrix.hpp"

namespace Ipopt
{

  SumSymMatrix::SumSymMatrix(const SumSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      factors_(owner_space->NTerms(), 1.0),
      matrices_(owner_space->NTerms()),
      owner_space_(owner_space)
  {}

  SumSymMatrix::~SumSymMatrix()
  {}

  void SumSymMatrix::SetTerm(Index iterm, Number factor,
                             const SymMatrix& matrix)
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factors_[iterm] = factor;
    matrices_[iterm] = &matrix;
  }

  void SumSymMatrix::GetTerm(Index iterm, Number& factor, SmartPtr<const SymMatrix>& matrix) const
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factor = factors_[iterm];
    matrix = matrices_[iterm];
  }


  Index SumSymMatrix::NTerms() const
  {
    return owner_space_->NTerms();
  }

  void SumSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                    Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      matrices_[iterm]->MultVector(alpha*factors_[iterm], x,
                                   1.0, y);
    }
  }

  void SumSymMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sSumSymMatrix \"%s\" of dimension %d with %d terms:\n",
            prefix.c_str(), name.c_str(), Dim(), NTerms());
    for (Index iterm=0; iterm<NTerms(); iterm++) {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "%sTerm %d with factor %23.16e and the following matrix:\n",
              prefix.c_str(), iterm, factors_[iterm]);
      char buffer[256];
      sprintf(buffer, "Term: %d", iterm);
      std::string name = buffer;
      matrices_[iterm]->Print(fp,name,indent,prefix);
    }
  }

  SumSymMatrix* SumSymMatrixSpace::MakeNewSumSymMatrix() const
  {
    return new SumSymMatrix(this);
  }

  SymMatrix* SumSymMatrixSpace::MakeNewSymMatrix() const
  {
    return MakeNewSumSymMatrix();
  }
} // namespace Ipopt
