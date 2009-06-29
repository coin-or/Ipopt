// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-05-12

#include "IpoptConfig.h"
#include "IpParVector.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

// FIXME - proper header files
//extern "C" {
#define MPICH_SKIP_MPICXX
#include "mpi.h"
//}

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParVector::ParVector(const ParVectorSpace* owner_space)
      :
      Vector(owner_space),
      owner_space_(owner_space)
  {
    DBG_START_METH("ParVector::ParVector()", dbg_verbosity);

    local_vector_ = owner_space_->LocalSpace()->MakeNewDenseVector();
  }

  ParVector::~ParVector()
  {
    DBG_START_METH("ParVector::~ParVector()", dbg_verbosity);
  }

  void ParVector::CopyImpl(const Vector& x)
  {
    DBG_START_METH("ParVector::CopyImpl(const Vector& x)", dbg_verbosity);
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->Copy(*p_x->LocalVector());
  }

  void ParVector::ScalImpl(Number alpha)
  {
    local_vector_->Scal(alpha);
  }

  void ParVector::AxpyImpl(Number alpha, const Vector &x)
  {
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->Axpy(alpha, *p_x->LocalVector());
  }

  Number ParVector::DotImpl(const Vector &x) const
  {
    Number retValue;
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    Number local_dot = local_vector_->Dot(*p_x->LocalVector());
    MPI_Allreduce(&local_dot, &retValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return retValue;
  }

  Number ParVector::Nrm2Impl() const
  {
    Number retValue;
    Number local_nrm2 = local_vector_->Nrm2();
    local_nrm2 *= local_nrm2;

    MPI_Allreduce(&local_nrm2, &retValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(retValue);
  }

  Number ParVector::AsumImpl() const
  {
    Number retValue;
    Number local_asum = local_vector_->Asum();

    MPI_Allreduce(&local_asum, &retValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return retValue;
  }

  Number ParVector::AmaxImpl() const
  {
    Number retValue;
    Number local_amax = local_vector_->Amax();

    MPI_Allreduce(&local_amax, &retValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return retValue;
  }

  void ParVector::SetImpl(Number value)
  {
    local_vector_->Set(value);
  }

  void ParVector::ElementWiseDivideImpl(const Vector& x)
  {
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->ElementWiseDivide(*p_x->LocalVector());
  }

  void ParVector::ElementWiseMultiplyImpl(const Vector& x)
  {
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->ElementWiseMultiply(*p_x->LocalVector());
  }

  void ParVector::ElementWiseMaxImpl(const Vector& x)
  {
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->ElementWiseMax(*p_x->LocalVector());
  }

  void ParVector::ElementWiseMinImpl(const Vector& x)
  {
    const ParVector* p_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    local_vector_->ElementWiseMin(*p_x->LocalVector());
  }

  void ParVector::ElementWiseReciprocalImpl()
  {
    local_vector_->ElementWiseReciprocal();
  }

  void ParVector::ElementWiseAbsImpl()
  {
    local_vector_->ElementWiseAbs();
  }

  void ParVector::ElementWiseSqrtImpl()
  {
    local_vector_->ElementWiseSqrt();
  }

  void ParVector::AddScalarImpl(Number scalar)
  {
    local_vector_->AddScalar(scalar);
  }

  Number ParVector::MaxImpl() const
  {
    Number retValue;
    Number local_max = local_vector_->Max();

    MPI_Allreduce(&local_max, &retValue, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return retValue;
  }

  Number ParVector::MinImpl() const
  {
    Number retValue;
    Number local_min = local_vector_->Min();

    MPI_Allreduce(&local_min, &retValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return retValue;
  }

  Number ParVector::SumImpl() const
  {
    Number retValue;
    Number local_sum = local_vector_->Sum();

    MPI_Allreduce(&local_sum, &retValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return retValue;
  }

  Number ParVector::SumLogsImpl() const
  {
    Number retValue;
    Number local_lsum = local_vector_->SumLogs();

    MPI_Allreduce(&local_lsum, &retValue, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return retValue;
  }

  void ParVector::ElementWiseSgnImpl()
  {
    local_vector_->ElementWiseSgn();
  }

  // Specialized Functions
  void ParVector::AddTwoVectorsImpl(Number a, const Vector& v1,
				    Number b, const Vector& v2, Number c)
  {
    const ParVector* p_v1 = static_cast<const ParVector*>(&v1);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&v1));

    const ParVector* p_v2 = static_cast<const ParVector*>(&v2);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&v2));

    local_vector_->AddTwoVectors(a, *p_v1->LocalVector(), b, *p_v2->LocalVector(), c);
  }

  Number
  ParVector::FracToBoundImpl(const Vector& delta, Number tau) const
  {
    Number retValue;
    const ParVector* p_delta = static_cast<const ParVector*>(&delta);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&delta));

    Number l_alpha = local_vector_->FracToBound(*p_delta->LocalVector(), tau);
    MPI_Allreduce(&l_alpha, &retValue, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    return retValue;
  }

  void ParVector::AddVectorQuotientImpl(Number a, const Vector& z,
                                          const Vector& s, Number c)
  {
    const ParVector* p_z = static_cast<const ParVector*>(&z);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&z));

    const ParVector* p_s = static_cast<const ParVector*>(&s);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&s));

    local_vector_->AddVectorQuotient(a, *p_z->LocalVector(), *p_s->LocalVector(), c);
  }

  void ParVector::PrintImpl(const Journalist& jnlst,
			    EJournalLevel level,
			    EJournalCategory category,
			    const std::string& name,
			    Index indent,
			    const std::string& prefix) const
  {
    if (Rank() == 0){
      jnlst.PrintfIndented(level, category, indent,
			   "%sParVector \"%s\" with %d pieces, dim %d:\n",
			   prefix.c_str(), name.c_str(), NumProc(), Dim());
    }
    char buffer[256];
    snprintf (buffer, 255, "%s[%d]", name.c_str(), Rank());
    std::string myname = buffer;
    
    local_vector_->Print( jnlst, level, category, myname, indent+1, prefix);
  }

  SmartPtr<DenseVector> ParVector::GlobalVectorNonConst() const
  {
    SmartPtr<DenseVector> globalv = owner_space_->GlobalSpace()->MakeNewDenseVector();

    Number* global_vals = globalv->Values();
    const Number* local_vals = local_vector_->Values();

    MPI_Allgatherv(const_cast<double*>(local_vals), LocalSize(),
		   MPI_DOUBLE, global_vals,
		   const_cast<int*>(&(owner_space_->RecvCounts()[0])),
		   const_cast<int*>(&(owner_space_->Displs()[0])),
		   MPI_DOUBLE, MPI_COMM_WORLD);

    return globalv;
  }

  SmartPtr<const DenseVector> ParVector::GlobalVector() const
  {
    SmartPtr<DenseVector> retval = GlobalVectorNonConst();
    return ConstPtr(retval);
  }  

  void ParVector::ExtractLocalVector(const DenseVector& global_vector)
  {
    local_vector_->CopyFromPos(owner_space_->StartPos(), global_vector);
    ObjectChanged();
  }

  void ParVector::ExtractLocalValues(const Number* vals)
  {
    local_vector_->SetValues(vals+owner_space_->StartPos());
    ObjectChanged();
  }

  ParVectorSpace::ParVectorSpace(Index total_dim, Index start_pos, Index size)
    :
    VectorSpace(total_dim),
    start_pos_(start_pos),
    size_(size)
  {
    local_space_ = new DenseVectorSpace(size);
    global_space_ = new DenseVectorSpace(total_dim);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc_);

    recvcounts_.resize(num_proc_);
    displs_.resize(num_proc_);

    MPI_Allgather(&start_pos_, 1, MPI_INT, &(displs_[0]), 1, MPI_INT, MPI_COMM_WORLD); 
    MPI_Allgather(&size_, 1, MPI_INT, &(recvcounts_[0]), 1, MPI_INT, MPI_COMM_WORLD); 

    // Test for consistency
    for (int i=0; i<num_proc_; i++){
      int nx = (i < num_proc_-1)? displs_[i+1] : total_dim;

      DBG_ASSERT(displs_[i] + recvcounts_[i] == nx);
    }

  }

} // namespace Ipopt
