// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpSymTMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"

namespace Ipopt
{

  SymTMatrix::SymTMatrix(const SymTMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space),
      values_(NULL),
      initialized_(false)
  {
    values_ = owner_space_->AllocateInternalStorage();

    if (Nonzeros() == 0) {
      initialized_ = true; // I guess ?!? what does this mean ?!?
    }
  }

  SymTMatrix::~SymTMatrix()
  {

    owner_space_->FreeInternalStorage(values_);
  }

  void SymTMatrix::SetValues(const Number* Values)
  {
    IpBlasDcopy(Nonzeros(), Values, 1, values_, 1);
    initialized_ = true;
    ObjectChanged();
  }

  void SymTMatrix::MultVectorImpl(Number alpha, const Vector &x, Number beta,
                                  Vector &y) const
  {
    //  A few sanity checks
    DBG_ASSERT(Dim()==x.Dim());
    DBG_ASSERT(Dim()==y.Dim());

    // Take care of the y part of the addition
    DBG_ASSERT(initialized_);
    if( beta!=0.0 ) {
      y.Scal(beta);
    }
    else {
      y.Set(0.0);  // In case y hasn't been initialized yet
    }

    // See if we can understand the data
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    DenseVector* dense_y = dynamic_cast<DenseVector*>(&y);
    DBG_ASSERT(dense_y); /* ToDo: Implement others */

    if (dense_x && dense_y) {
      const Index*  irn=Irows();
      const Index*  jcn=Jcols();
      const Number* val=values_;
      const Number* xvals=dense_x->Values();
      Number* yvals=dense_y->Values();
      for(Index i=0; i<Nonzeros(); i++) {
        yvals[*irn-1] += alpha* (*val) * xvals[*jcn-1];
        if (*irn!=*jcn) {
          // this is not a diagonal element
          yvals[*jcn-1] += alpha* (*val) * xvals[*irn-1];
        }
        val++;
        irn++;
        jcn++;
      }
    }
  }

  Number* SymTMatrix::Values()
  {
    // cannot check for initialized values here, in case this pointer is
    // requested for setting the first values
    //DBG_ASSERT(initialized_);

    // Here we assume that every time someone requests this direct raw
    // pointer, the data is going to change and the Tag for this
    // vector has to be updated.
    ObjectChanged();
    initialized_ = true;
    return values_;
  }

  const Number* SymTMatrix::Values() const
  {
    DBG_ASSERT(initialized_);
    return values_;
  }

  void SymTMatrix::FillStruct(ipfint* Irn, ipfint* Jcn) const
  {
    DBG_ASSERT(initialized_);
    for(Index i=0; i<Nonzeros(); i++) {
      Irn[i] = Irows()[i];
      Jcn[i] = Jcols()[i];
    }
  }

  void SymTMatrix::FillValues(Number* Values) const
  {
    DBG_ASSERT(initialized_);
    IpBlasDcopy(Nonzeros(), values_, 1, Values, 1);
  }

  void SymTMatrix::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    fprintf(fp, "\n");
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sSymTMatrix \"%s\" with %d nonzero elements:\n",
            prefix.c_str(), name.c_str(), Nonzeros());
    if (initialized_) {
      for (Index i=0; i<Nonzeros(); i++) {
        for (Index ind=0; ind<indent; ind++) {
          fprintf(fp, " ");
        }
        fprintf(fp, "%s%s[%5d,%5d]=%23.16e  (%d)\n", prefix.c_str(), name.c_str(), Irows()[i],
                Jcols()[i], values_[i], i);
      }
    }
    else {
      for (Index ind=0; ind<indent; ind++) {
        fprintf(fp, " ");
      }
      fprintf(fp, "%sUninitialized! ", prefix.c_str());
    }
  }

  SymTMatrixSpace::SymTMatrixSpace(Index dim, Index nonZeros,
                                   const Index* iRows,
                                   const Index* jCols)
      :
      SymMatrixSpace(dim),
      nonZeros_(nonZeros),
      iRows_(NULL),
      jCols_(NULL)
  {
    iRows_ = new Index[nonZeros];
    jCols_ = new Index[nonZeros];
    for (Index i=0; i<nonZeros; i++) {
      iRows_[i] = iRows[i];
      jCols_[i] = jCols[i];
    }
  }

  SymTMatrixSpace::~SymTMatrixSpace()
  {
    delete [] iRows_;
    delete [] jCols_;
  }

  Number* SymTMatrixSpace::AllocateInternalStorage() const
  {
    return new Number[Nonzeros()];
  }

  void SymTMatrixSpace::FreeInternalStorage(Number* values) const
  {
    delete [] values;
  }

} // namespace Ipopt
