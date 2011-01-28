// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpGenTMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace Ipopt
{

  GenTMatrix::GenTMatrix(const GenTMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space),
      values_(NULL),
      initialized_(false)
  {
    values_ = owner_space_->AllocateInternalStorage();

    if (Nonzeros() == 0) {
      initialized_ = true; // I guess ?!? what does this mean ?!?
    }
  }

  GenTMatrix::~GenTMatrix()
  {
    owner_space_->FreeInternalStorage(values_);
  }

  void GenTMatrix::SetValues(const Number* Values)
  {
    IpBlasDcopy(Nonzeros(), Values, 1, values_, 1);
    initialized_ = true;
    ObjectChanged();
  }

  void GenTMatrix::MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                  Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==x.Dim());
    DBG_ASSERT(NRows()==y.Dim());

    // Take care of the y part of the addition
    DBG_ASSERT(initialized_);
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
    DenseVector* dense_y = static_cast<DenseVector*>(&y);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

    if (dense_x && dense_y) {
      const Index*  irows=Irows();
      const Index*  jcols=Jcols();
      const Number* val=values_;
      Number* yvals=dense_y->Values();
      yvals--;
      if (dense_x->IsHomogeneous()) {
        Number as = alpha * dense_x->Scalar();
        for (Index i=0; i<Nonzeros(); i++) {
          yvals[*irows] += as * (*val);
          val++;
          irows++;
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        xvals--;
        for (Index i=0; i<Nonzeros(); i++) {
          yvals[*irows] += alpha* (*val) * xvals[*jcols];
          val++;
          irows++;
          jcols++;
        }
      }
    }
  }

  void GenTMatrix::TransMultVectorImpl(Number alpha, const Vector &x, Number beta,
                                       Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(NCols()==y.Dim());
    DBG_ASSERT(NRows()==x.Dim());

    // Take care of the y part of the addition
    DBG_ASSERT(initialized_);
    if ( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = static_cast<const DenseVector*>(&x);
    DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));
    DenseVector* dense_y = static_cast<DenseVector*>(&y);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&y));

    if (dense_x && dense_y) {
      const Index*  irows=Irows();
      const Index*  jcols=Jcols();
      const Number* val=values_;
      Number* yvals=dense_y->Values();
      yvals--;

      if (dense_x->IsHomogeneous()) {
        Number as = alpha * dense_x->Scalar();
        for (Index i=0; i<Nonzeros(); i++) {
          yvals[*jcols] += as * (*val);
          val++;
          jcols++;
        }
      }
      else {
        const Number* xvals=dense_x->Values();
        xvals--;
        for (Index i=0; i<Nonzeros(); i++) {
          yvals[*jcols] += alpha* (*val) * xvals[*irows];
          val++;
          irows++;
          jcols++;
        }
      }
    }
  }

  bool GenTMatrix::HasValidNumbersImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum = IpBlasDasum(Nonzeros(), values_, 1);
    return IsFiniteNumber(sum);
  }

  void GenTMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    DBG_ASSERT(initialized_);

    DenseVector* dense_vec = static_cast<DenseVector*>(&rows_norms);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&rows_norms));

    const Index* irows=Irows();
    const Number* val=values_;
    Number* vec_vals=dense_vec->Values();
    vec_vals--;

    for (Index i=0; i<Nonzeros(); i++) {
      vec_vals[irows[i]] = Max(vec_vals[irows[i]], fabs(val[i]));
    }
  }

  void GenTMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
    DBG_ASSERT(initialized_);

    DenseVector* dense_vec = static_cast<DenseVector*>(&cols_norms);
    DBG_ASSERT(dynamic_cast<DenseVector*>(&cols_norms));

    const Index* jcols=Jcols();
    const Number* val=values_;
    Number* vec_vals=dense_vec->Values();
    vec_vals--;

    for (Index i=0; i<Nonzeros(); i++) {
      vec_vals[jcols[i]] = Max(vec_vals[jcols[i]], fabs(val[i]));
    }
  }

  void GenTMatrix::PrintImplOffset(const Journalist& jnlst,
                                   EJournalLevel level,
                                   EJournalCategory category,
                                   const std::string& name,
                                   Index indent,
                                   const std::string& prefix,
                                   Index offset) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sGenTMatrix \"%s\" of dimension %d by %d with %d nonzero elements:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols(), Nonzeros());
    if (initialized_) {
      for (Index i=0; i<Nonzeros(); i++) {
        jnlst.PrintfIndented(level, category, indent,
                             "%s%s[%5d,%5d]=%23.16e  (%d)\n",
                             prefix.c_str(), name.c_str(), Irows()[i]+offset,
                             Jcols()[i], values_[i], i);
      }
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sUninitialized!\n", prefix.c_str());
    }
  }

  GenTMatrixSpace::GenTMatrixSpace(Index nRows, Index nCols,
                                   Index nonZeros,
                                   const Index* iRows, const Index* jCols)
      :
      MatrixSpace(nRows, nCols),
      nonZeros_(nonZeros),
      jCols_(NULL),
      iRows_(NULL)
  {
    iRows_ = new Index[nonZeros];
    jCols_ = new Index[nonZeros];
    for (Index i=0; i<nonZeros; i++) {
      iRows_[i] = iRows[i];
      jCols_[i] = jCols[i];
    }
  }

  Number* GenTMatrixSpace::AllocateInternalStorage() const
  {
    return new Number[Nonzeros()];
  }

  void GenTMatrixSpace::FreeInternalStorage(Number* values) const
  {
    delete [] values;
  }


} // namespace Ipopt
