// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMatrix.hpp"
#include "IpIdentityMatrix.hpp"
#include "IpZeroMatrix.hpp"

namespace Ipopt
{

  Matrix::Matrix(const MatrixSpace* owner_space)
      :
      TaggedObject(),
      owner_space_(owner_space)
  {}

  Matrix::~Matrix()
  {}

  void Matrix::MultVector(Number alpha, const Vector& x, Number beta,
                          Vector& y) const
  {
    MultVectorImpl(alpha, x, beta, y);
  }

  void Matrix::TransMultVector(Number alpha, const Vector& x, Number beta,
                               Vector& y) const
  {
    TransMultVectorImpl(alpha, x, beta, y);
  }

  void Matrix::AddMSinvZ(Number alpha, const Vector& S, const Vector& Z,
                         Vector& X) const
  {
    AddMSinvZImpl(alpha, S, Z, X);
  }

  void Matrix::SinvBlrmZMTdBr(Number alpha, const Vector& S,
                              const Vector& R, const Vector& Z,
                              const Vector& D, Vector& X) const
  {
    SinvBlrmZMTdBrImpl(alpha, S, R, Z, D, X);
  }

  // Prototype for specialize methods (can and should be overloaded)
  void Matrix::AddMSinvZImpl(Number alpha, const Vector& S, const Vector& Z,
                             Vector& X) const
  {
    SmartPtr<Vector> tmp = S.MakeNew();
    tmp->AddVectorQuotient(1., Z, S, 0.);
    MultVector(alpha, *tmp, 1., X);
  }

  void Matrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
                                  const Vector& R, const Vector& Z,
                                  const Vector& D, Vector& X) const
  {
    TransMultVector(alpha, D, 0., X);
    X.ElementWiseMultiply(Z);
    X.Axpy(1., R);
    X.ElementWiseDivide(S);
  }

  void Matrix::Print(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    PrintImpl(fp, name, indent, prefix);
  }

  MatrixSpace::MatrixSpace(Index nRows, Index nCols)
      :
      nRows_(nRows),
      nCols_(nCols)
  {}

  MatrixSpace::~MatrixSpace()
  {}

  bool MatrixSpace::IsMatrixFromSpace(const Matrix& matrix) const
  {
    return (matrix.OwnerSpace() == this);
  }

} // namespace Ipopt
