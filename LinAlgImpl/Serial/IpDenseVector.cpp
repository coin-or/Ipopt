// Copyright (C) 2004, International BusinDess Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpDenseVector.hpp"
#include "IpBlas.hpp"
#include "IpUtils.hpp"
#include "IpDebug.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

namespace Ipopt
{

  //  static const Index dbg_verbosity = 0;

  DenseVector::DenseVector(const DenseVectorSpace* owner_space)
      :
      Vector(owner_space),
      owner_space_(owner_space),
      values_(NULL),
      initialized_(false)
  {
    DBG_START_METH("DenseVector::DenseVector(Index dim)", dbg_verbosity);
    values_ = owner_space_->AllocateInternalStorage();
    if (Dim() == 0) {
      initialized_ = true; // I guess ?!? what does this even mean ?!?
    }
  }

  DenseVector::~DenseVector()
  {
    DBG_START_METH("DenseVector::~DenseVector()", dbg_verbosity);
    owner_space_->FreeInternalStorage(values_);
  }

  void DenseVector::SetValues(const Number* x)
  {
    initialized_ = true;
    IpBlasDcopy(Dim(), x, 1, values_, 1);
    // This is not an overloaded method from
    // Vector. Here, we must call ObjectChanged()
    // manually.
    ObjectChanged();
  }

  Number* DenseVector::Values()
  {
    // Here we assume that every time someone requests this direct raw
    // pointer, the data is going to change and the Tag for this
    // vector has to be updated.
    ObjectChanged();
    initialized_= true;
    return values_;
  }

  const Number* DenseVector::Values() const
  {
    DBG_ASSERT(initialized_);
    return values_;
  }

  void DenseVector::CopyImpl(const Vector& x)
  {
    DBG_START_METH("DenseVector::CopyImpl(const Vector& x)", dbg_verbosity);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x); /* ToDo: Implement others */
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      IpBlasDcopy(Dim(), dense_x->values_, 1, values_, 1);
    }
    initialized_=true;
  }

  void DenseVector::ScalImpl(Number alpha)
  {
    DBG_ASSERT(initialized_);
    IpBlasDscal(Dim(), alpha, values_, 1);
  }

  void DenseVector::AxpyImpl(Number alpha, const Vector &x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      IpBlasDaxpy(Dim(), alpha, dense_x->values_, 1, values_, 1);
    }
  }

  Number DenseVector::DotImpl(const Vector &x) const
  {
    DBG_ASSERT(initialized_);
    Number retValue = 0.0;
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      retValue = IpBlasDdot(Dim(), dense_x->values_, 1, values_, 1);
    }
    return retValue;
  }

  Number DenseVector::Nrm2Impl() const
  {
    DBG_ASSERT(initialized_);
    return IpBlasDnrm2(Dim(), values_, 1);
  }

  Number DenseVector::AsumImpl() const
  {
    DBG_ASSERT(initialized_);
    return IpBlasDasum(Dim(), values_, 1);
  }

  Number DenseVector::AmaxImpl() const
  {
    DBG_ASSERT(initialized_);
    if (Dim()==0) {
      return 0.;
    }
    else {
      return fabs(values_[IpBlasIdamax(Dim(), values_, 1)-1]);
    }
  }

  void DenseVector::SetImpl(Number value)
  {
    initialized_ = true;
    IpBlasDcopy(Dim(), &value, 0, values_, 1);
  }

  void DenseVector::ElementWiseDivideImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      // Is there a BLAS type operation for an efficient implementation? NO
      for (Index i=0; i<Dim(); i++) {
        values_[i] = values_[i]/dense_x->values_[i];
      }
    }
  }

  void DenseVector::ElementWiseMultiplyImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      // Is there a BLAS type operation for an efficient implementation? NO
      for (Index i=0; i<Dim(); i++) {
        values_[i] = values_[i]*dense_x->values_[i];
      }
    }
  }

  void DenseVector::ElementWiseMaxImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      for (Index i=0; i<Dim(); i++) {
        values_[i] = Ipopt::Max(values_[i], dense_x->values_[i]);
      }
    }
  }

  void DenseVector::ElementWiseMinImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(Dim() == dense_x->Dim());
      for (Index i=0; i<Dim(); i++) {
        values_[i] = Ipopt::Min(values_[i], dense_x->values_[i]);
      }
    }
  }

  void DenseVector::ElementWiseReciprocalImpl()
  {
    DBG_ASSERT(initialized_);
    for (Index i=0; i<Dim(); i++) {
      values_[i] = 1.0/values_[i];
    }
  }

  void DenseVector::ElementWiseSqrtImpl()
  {
    DBG_ASSERT(initialized_);
    for (Index i=0; i<Dim(); i++) {
      values_[i] = sqrt(values_[i]);
    }
  }

  void DenseVector::AddScalarImpl(Number scalar)
  {
    DBG_ASSERT(initialized_);
    IpBlasDaxpy(Dim(), 1., &scalar, 0, values_, 1);
  }

  Number DenseVector::MaxImpl() const
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(Dim() > 0 && "There is no Max of a zero length vector (no reasonable default can be returned)");
    Number max = 0.0;
    if (Dim() > 0) {
      max = values_[0];
      for (Index i=1; i<Dim(); i++) {
        max = Ipopt::Max(values_[i], max);
      }
    }
    return max;
  }

  Number DenseVector::MinImpl() const
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(Dim() > 0 && "There is no Min of a zero length vector"
               "(no reasonable default can be returned) - "
               "Check for zero length vector before calling");
    Number min = 0.0;
    if (Dim() > 0) {
      min = values_[0];
      for (Index i=1; i<Dim(); i++) {
        min = Ipopt::Min(values_[i], min);
      }
    }
    return min;
  }

  Number DenseVector::SumImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum = 0.0;
    for (Index i=0; i<Dim(); i++) {
      sum += values_[i];
    }
    return sum;
  }

  Number DenseVector::SumLogsImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum = 0.0;
    for (Index i=0; i<Dim(); i++) {
      sum += log(values_[i]);
    }
    return sum;
  }

  void DenseVector::SgnImpl()
  {
    DBG_ASSERT(initialized_);
    for (Index i=0; i<Dim(); i++) {
      if (values_[i] > 0) {
        values_[i] = 1;
      }
      else if (values_[i] < 0) {
        values_[i] = -1;
      }
      else {
        values_[i] = 0;
      }
    }
  }

  void DenseVector::CopyToPos(Index Pos, const Vector& x)
  {
    Index dim_x = x.Dim();
    DBG_ASSERT(dim_x+Pos<=Dim());
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    IpBlasDcopy(dim_x, dense_x->values_, 1, values_+Pos, 1);
    initialized_=true;
    ObjectChanged();
  }

  void DenseVector::CopyFromPos(Index Pos, Vector& x) const
  {
    Index dim_x = x.Dim();
    DBG_ASSERT(dim_x+Pos<=Dim());
    DenseVector* dense_x = dynamic_cast<DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    IpBlasDcopy(dim_x, values_+Pos, 1, dense_x->values_, 1);
    // We need to tell X that it has changed!
    dense_x->ObjectChanged();
    dense_x->initialized_=true;
  }

  void DenseVector::PrintImpl(FILE* fp, std::string name, Index indent, std::string prefix) const
  {
    for (Index ind=0; ind<indent; ind++) {
      fprintf(fp, " ");
    }
    fprintf(fp, "%sDenseVector \"%s\" with %d elements:\n",
            prefix.c_str(), name.c_str(), Dim());
    if (initialized_) {
      for (Index i=0; i<Dim(); i++) {
        for (Index ind=0; ind<indent; ind++) {
          fprintf(fp, " ");
        }
        fprintf(fp, "%s%s[%5d]=%23.16e\n", prefix.c_str(), name.c_str(), i+1, values_[i]);
      }
    }
    else {
      fprintf(fp, "Uninitialized!\n");
    }
  }

  Number* DenseVectorSpace::AllocateInternalStorage() const
  {
    return new Number[Dim()];
  }

  void DenseVectorSpace::FreeInternalStorage(Number* values) const
  {
    delete [] values;
  }

} // namespace Ipopt
