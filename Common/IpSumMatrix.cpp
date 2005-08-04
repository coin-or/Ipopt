// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSumMatrix.hpp"

namespace Ipopt
{

  SumMatrix::SumMatrix(const SumMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      factors_(owner_space->NTerms(), 1.0),
      matrices_(owner_space->NTerms()),
      owner_space_(owner_space)
  {}

  SumMatrix::~SumMatrix()
  {}

  void SumMatrix::SetTerm(Index iterm, Number factor,
                          const Matrix& matrix)
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factors_[iterm] = factor;
    matrices_[iterm] = &matrix;
  }

  void SumMatrix::GetTerm(Index iterm, Number& factor, SmartPtr<const Matrix>& matrix) const
  {
    DBG_ASSERT(iterm<owner_space_->NTerms());
    factor = factors_[iterm];
    matrix = matrices_[iterm];
  }

  Index SumMatrix::NTerms() const
  {
    return owner_space_->NTerms();
  }

  void SumMatrix::MultVectorImpl(Number alpha, const Vector &x,
                                 Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());

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

  void SumMatrix::TransMultVectorImpl(Number alpha, const Vector& x,
                                      Number beta, Vector& y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NRows()==x.Dim());
    DBG_ASSERT(NCols()==y.Dim());

    // Take care of the y part of the addition
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    for (Index iterm=0; iterm<NTerms(); iterm++) {
      DBG_ASSERT(IsValid(matrices_[iterm]));
      matrices_[iterm]->TransMultVector(alpha*factors_[iterm], x,
                                        1.0, y);
    }
  }



  void SumMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sSumMatrix \"%s\" of dimension %d x %d with %d terms:\n",
            prefix.c_str(), name.c_str(), NRows(), NCols(), NTerms());
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

  SumMatrix* SumMatrixSpace::MakeNewSumMatrix() const
  {
    return new SumMatrix(this);
  }

  Matrix* SumMatrixSpace::MakeNew() const
  {
    return MakeNewSumMatrix();
  }
} // namespace Ipopt
