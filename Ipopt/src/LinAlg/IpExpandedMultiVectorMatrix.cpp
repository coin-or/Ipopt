// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter             IBM    2009-11-05

#include "IpExpandedMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"

#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ExpandedMultiVectorMatrix::ExpandedMultiVectorMatrix(const ExpandedMultiVectorMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space),
      vecs_(owner_space->NRows())
  {}

  void ExpandedMultiVectorMatrix::SetVector(Index i, SmartPtr<const Vector> vec)
  {
    DBG_ASSERT(i<NRows());
    vecs_[i] = vec;
    ObjectChanged();
  }

  void ExpandedMultiVectorMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NRows()==x.Dim());
    DBG_ASSERT(NCols()==y.Dim());

    SmartPtr<const ExpansionMatrix> P = GetExpansionMatrix();
    SmartPtr<Vector> y_tmp;

    if (IsValid(P)) {
      y_tmp = RowVectorSpace()->MakeNew();
      y_tmp->Set(0.);
    }
    else {
      // Take care of the y part of the addition
      if ( beta!=0.0 ) {
        y.Scal(beta);
      }
      else {
        y.Set(0.0);  // In case y hasn't been initialized yet
      }
      y_tmp = &y;
    }

    // See if we can understand the data
    const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));

    // We simply add all the Vectors one after the other
    if (dense_x->IsHomogeneous()) {
      Number val = dense_x->Scalar();
      for (Index i=0; i<NRows(); i++) {
        if (IsValid(vecs_[i])) {
          y_tmp->AddOneVector(alpha*val, *vecs_[i], 1.);
        }
      }
    }
    else {
      const Number* values = dense_x->Values();
      for (Index i=0; i<NRows(); i++) {
        if (IsValid(vecs_[i])) {
          y_tmp->AddOneVector(alpha*values[i], *vecs_[i], 1.);
        }
      }
    }

    if (IsValid(P)) {
      P->MultVector(1., *y_tmp, beta, y);
    }
  }

  void ExpandedMultiVectorMatrix::MultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NRows()==y.Dim());
    DBG_ASSERT(NCols()==x.Dim());

    // check if there is an expansion matrix
    SmartPtr<const ExpansionMatrix> P = GetExpansionMatrix();
    SmartPtr<const Vector> x_tmp;
    if (IsValid(P)) {
      SmartPtr<Vector> exp_x = RowVectorSpace()->MakeNew();
      P->TransMultVector(1., x, 0., *exp_x);
      x_tmp = ConstPtr(exp_x);
    }
    else {
      x_tmp = &x;
    }

    // See if we can understand the data
    DenseVector* dense_y = static_cast<DenseVector*>(&y);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

    // Use the individual dot products to get the matrix (transpose)
    // vector product
    Number *yvals=dense_y->Values();
    if ( beta!=0.0 ) {
      for (Index i=0; i<NRows(); i++) {
        if (IsValid(vecs_[i])) {
          yvals[i] = alpha*vecs_[i]->Dot(*x_tmp) + beta*yvals[i];
        }
        else {
          yvals[i] *= beta;
        }
      }
    }
    else {
      for (Index i=0; i<NRows(); i++) {
        if (IsValid(vecs_[i])) {
          yvals[i] = alpha*vecs_[i]->Dot(*x_tmp);
        }
        else {
          yvals[i] = 0.;
        }
      }
    }
  }

  bool ExpandedMultiVectorMatrix::HasValidNumbersImpl() const
  {
    for (Index i=0; i<NRows(); i++) {
      if (IsValid(vecs_[i]) && vecs_[i]->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void
  ExpandedMultiVectorMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    THROW_EXCEPTION(UNIMPLEMENTED_LINALG_METHOD_CALLED,
                    "ExpandedMultiVectorMatrix::ComputeRowAMaxImpl not implemented");
  }

  void
  ExpandedMultiVectorMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
    THROW_EXCEPTION(UNIMPLEMENTED_LINALG_METHOD_CALLED,
                    "ExpandedMultiVectorMatrix::ComputeColAMaxImpl not implemented");
  }

  void ExpandedMultiVectorMatrix::PrintImpl(const Journalist& jnlst,
      EJournalLevel level,
      EJournalCategory category,
      const std::string& name,
      Index indent,
      const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sExpandedMultiVectorMatrix \"%s\" with %d columns:\n",
                         prefix.c_str(), name.c_str(), NRows());

    for (Index i=0; i<NRows(); i++) {
      if (IsValid(vecs_[i])) {
        DBG_ASSERT(name.size()<200);
        char buffer[256];
        Snprintf(buffer, 255, "%s[%2d]", name.c_str(), i);
        std::string term_name = buffer;
        vecs_[i]->Print(&jnlst, level, category, term_name,
                        indent+1, prefix);
      }
      else {
        jnlst.PrintfIndented(level, category, indent,
                             "%sVector in column %d is not yet set!\n",
                             prefix.c_str(), i);
      }
    }
    SmartPtr<const ExpansionMatrix> P = GetExpansionMatrix();
    if (IsValid(P)) {
      char buffer[256];
      Snprintf(buffer, 255, "%s[ExpMat]", name.c_str());
      std::string term_name = buffer;
      P->Print(&jnlst, level, category, term_name,
               indent+1, prefix);
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sExpandedMultiVectorMatrix \"%s\" has no ExpansionMatrix\n", prefix.c_str(), name.c_str());
    }
  }

  ExpandedMultiVectorMatrixSpace::
  ExpandedMultiVectorMatrixSpace(Index nrows,
                                 const VectorSpace& vec_space,
                                 SmartPtr<const ExpansionMatrix> exp_matrix)
      :
      MatrixSpace(nrows, IsValid(exp_matrix) ? exp_matrix->NRows() : vec_space.Dim()),
      vec_space_(&vec_space),
      exp_matrix_(exp_matrix)
  {}

} // namespace Ipopt
