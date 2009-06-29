// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#include "IpParExpansionMatrix.hpp"
#include "IpParVector.hpp"

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

  ParExpansionMatrix::ParExpansionMatrix(const ParExpansionMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {
    DBG_START_METH("ParGenMatrix::ParGenMatrix()", dbg_verbosity);

    local_matrix_ = owner_space_->LocalSpace()->MakeNewExpansionMatrix();
  }

  ParExpansionMatrix::~ParExpansionMatrix()
  {
    DBG_START_METH("ParGenMatrix::~ParGenMatrix()", dbg_verbosity);
  }

  void ParExpansionMatrix::MultVectorImpl(Number alpha, const Vector &x,
					  Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));
    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

    local_matrix_->MultVector(alpha, *par_x->LocalVector(),
			      beta, *par_y->LocalVector());
  }

  void ParExpansionMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));
    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

    local_matrix_->TransMultVector(alpha, *par_x->LocalVector(),
				   beta, *par_y->LocalVector());
  }

  // Specialized method (overloaded from IpMatrix)
  void ParExpansionMatrix::AddMSinvZImpl(Number alpha, const Vector& S,
                                      const Vector& Z, Vector& X) const
  {
    const ParVector* par_S = static_cast<const ParVector*>(&S);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&S));
    const ParVector* par_Z = static_cast<const ParVector*>(&Z);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&Z));
    ParVector* par_X = static_cast<ParVector*>(&X);
    DBG_ASSERT(dynamic_cast<ParVector*>(&X));

    local_matrix_->AddMSinvZ(alpha, *par_S->LocalVector(),
			     *par_Z->LocalVector(), *par_X->LocalVector());
  }

  void ParExpansionMatrix::SinvBlrmZMTdBrImpl(Number alpha, const Vector& S,
					      const Vector& R, const Vector& Z,
					      const Vector& D, Vector& X) const
  {
    DBG_START_METH("ParExpansionMatrix::SinvBlrmZMTdBrImpl", dbg_verbosity);

    const ParVector* par_S = static_cast<const ParVector*>(&S);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&S));
    const ParVector* par_R = static_cast<const ParVector*>(&R);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&R));
    const ParVector* par_Z = static_cast<const ParVector*>(&Z);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&Z));
    const ParVector* par_D = static_cast<const ParVector*>(&D);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&D));
    ParVector* par_X = static_cast<ParVector*>(&X);
    DBG_ASSERT(dynamic_cast<ParVector*>(&X));

    local_matrix_->SinvBlrmZMTdBr(alpha, *par_S->LocalVector(),
				  *par_R->LocalVector(), *par_Z->LocalVector(),
				  *par_D->LocalVector(), *par_X->LocalVector());
  }

  void ParExpansionMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    ParVector* par_vec = static_cast<ParVector*>(&rows_norms);
    DBG_ASSERT(dynamic_cast<ParVector*>(&rows_norms));

    local_matrix_->ComputeRowAMax(*par_vec->LocalVector(), init);
  }

  void ParExpansionMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
    ParVector* par_vec = static_cast<ParVector*>(&cols_norms);
    DBG_ASSERT(dynamic_cast<ParVector*>(&cols_norms));

    local_matrix_->ComputeColAMax(*par_vec->LocalVector(), init);
  }

  void ParExpansionMatrix::PrintImpl(const Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category,
                                  const std::string& name,
                                  Index indent,
                                  const std::string& prefix) const
  {
    const int rank = owner_space_->SmallVectorSpace()->Rank();
    if (rank == 0){
      jnlst.PrintfIndented(level, category, indent,
			   "%sParExpansionMatrix \"%s\"\n",
			   prefix.c_str(), name.c_str());
    }

    char buffer[256];
    snprintf (buffer, 255, "%s[%d]", name.c_str(), rank);
    std::string myname = buffer;

    local_matrix_->Print(jnlst, level, category, myname, indent+1, prefix);
  }

  ParExpansionMatrixSpace::
  ParExpansionMatrixSpace(SmartPtr<const ParVectorSpace> LargeVectorSpace, 
			  SmartPtr<const ParVectorSpace> SmallVectorSpace, 
			  const Index *ExpPos,
			  const int offset)
    :
    MatrixSpace(LargeVectorSpace->LocalSize(), SmallVectorSpace->LocalSize()),
    large_vector_space_(LargeVectorSpace),
    small_vector_space_(SmallVectorSpace),
    global_compressed_pos_(NULL)
  {
    local_space_ = new ExpansionMatrixSpace(LargeVectorSpace->LocalSize(),
					    SmallVectorSpace->LocalSize(),
					    ExpPos, offset);
  }

  ParExpansionMatrixSpace::~ParExpansionMatrixSpace()
  {
    delete [] global_compressed_pos_;
  }

  const Index* ParExpansionMatrixSpace::
  GlobalCompressedPosIndices() const
  {
    if (!global_compressed_pos_) {
      global_compressed_pos_ = new Index[large_vector_space_->Dim()];

      const Index local_nrows = large_vector_space_->LocalSpace()->Dim();
      Index* local_compressed_pos = new Index[local_nrows];
      const Index* orig_pos = LocalSpace()->CompressedPosIndices();
      const Index offset = small_vector_space_->StartPos();
      for (Index i=0; i<local_nrows; i++) {
	if (orig_pos[i] == -1) {
	  local_compressed_pos[i] = -1;
	}
	else {
	  local_compressed_pos[i] = orig_pos[i] + offset;
	}
      }
      MPI_Allgatherv(local_compressed_pos, local_nrows,
		     MPI_INT, global_compressed_pos_,
		     const_cast<int*>(&(large_vector_space_->RecvCounts()[0])),
		     const_cast<int*>(&(large_vector_space_->Displs()[0])),
		     MPI_INT, MPI_COMM_WORLD);

      delete [] local_compressed_pos;
    }
    return global_compressed_pos_;
  }

} // namespace Ipopt
