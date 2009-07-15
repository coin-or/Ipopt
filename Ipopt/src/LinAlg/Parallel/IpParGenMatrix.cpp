// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-05-29

#include "IpParVector.hpp"
#include "IpParGenMatrix.hpp"
#include "IpBlas.hpp"
#include "IpLapack.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include "IpMpi.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParGenMatrix::ParGenMatrix(const ParGenMatrixSpace* owner_space)
      :
      Matrix(owner_space),
      owner_space_(owner_space)
  {
    DBG_START_METH("ParGenMatrix::ParGenMatrix()", dbg_verbosity);

    local_matrix_ = owner_space_->LocalSpace()->MakeNewGenTMatrix();
  }

  ParGenMatrix::~ParGenMatrix()
  {
    DBG_START_METH("ParGenMatrix::~ParGenMatrix()", dbg_verbosity);
  }

  // assume x is parallel, and y is parallel
  void ParGenMatrix::MultVectorImpl(Number alpha, const Vector &x,
				    Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

    SmartPtr<const DenseVector> dense_x = par_x->GlobalVector();
    SmartPtr<DenseVector> local_y = par_y->LocalVector();

    local_matrix_->MultVector(alpha, *dense_x, beta, *local_y);
  }

  // assume x is parallel, and y is parallel
  void ParGenMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
					 Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

    SmartPtr<const DenseVector> local_x = par_x->LocalVector();
    SmartPtr<DenseVector> dense_y = par_y->MakeNewGlobalVector();

    local_matrix_->TransMultVector(alpha, *local_x, 0., *dense_y);

    // CAN THIS be more efficient?
    Number *yvalues = dense_y->Values();
    MPI_Allreduce(MPI_IN_PLACE, yvalues, NCols(),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (beta==0.) {
      par_y->ExtractLocalVector(*dense_y);
    }
    else {
      SmartPtr<ParVector> par_v = par_y->MakeNewParVector();
      par_v->ExtractLocalVector(*dense_y);
      par_y->Axpy(beta, *par_v);
    }
  }

  void ParGenMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    ParVector* par_vec = static_cast<ParVector*>(&rows_norms);
    DBG_ASSERT(dynamic_cast<ParVector*>(&rows_norms));
    DenseVector* local_vec =par_vec->LocalVector();

    local_matrix_->ComputeRowAMax(*local_vec, init);
  }

  void ParGenMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
    DBG_ASSERT(dynamic_cast<DenseVector*>(&cols_norms));

    throw("ParGenMatrix::ComputeColAMaxImpl not properly implemented");
    // This will crash in the GenTMatrix, because cols_norms is not a
    // DenseVector
    local_matrix_->ComputeColAMax(cols_norms, init);

    DenseVector* dense_vec = static_cast<DenseVector*>(&cols_norms);
    Number* vec_vals=dense_vec->Values();

    MPI_Allreduce(MPI_IN_PLACE, vec_vals, NCols(),
		  MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }

  bool ParGenMatrix::HasValidNumbersImpl() const
  {
    int valid = local_matrix_->HasValidNumbers();
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    return (bool)valid;
  }

  void ParGenMatrix::PrintImpl(const Journalist& jnlst,
			       EJournalLevel level,
			       EJournalCategory category,
			       const std::string& name,
			       Index indent,
			       const std::string& prefix) const
  {
    jnlst.PrintfIndented(level, category, indent,
			 "%sParVector \"%s\" with %d pieces, nrows %d, ncols:\n",
			 prefix.c_str(), name.c_str(), NumProc(), NRows(), NCols());
    char buffer[256];
    snprintf (buffer, 255, "%s[%2d]", name.c_str(), Rank());
    std::string myname = buffer;
    
    jnlst.StartDistributedOutput();
    local_matrix_->PrintImplOffset(jnlst, level, category, myname, indent+1,
				   prefix, RowStartPos());
    jnlst.FinishDistributedOutput();
  }

  ParGenMatrixSpace::ParGenMatrixSpace(SmartPtr<const ParVectorSpace> RowVectorSpace, Index nCols,
				       Index nonZeros,
				       const Index* iRows, const Index* jCols)
    :
    MatrixSpace(RowVectorSpace->Dim(), nCols),
    rowVectorSpace_(RowVectorSpace)
  {
    local_space_ = new GenTMatrixSpace(rowVectorSpace_->LocalSize(),
				       nCols, nonZeros, iRows, jCols);
  }

} // namespace Ipopt
