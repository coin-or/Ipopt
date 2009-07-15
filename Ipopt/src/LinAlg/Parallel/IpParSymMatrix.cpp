// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-05-29

#include "IpParVector.hpp"
#include "IpParSymMatrix.hpp"
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

  ParSymMatrix::ParSymMatrix(const ParSymMatrixSpace* owner_space)
      :
      SymMatrix(owner_space),
      owner_space_(owner_space)
  {
    DBG_START_METH("ParSymMatrix::ParSymMatrix()", dbg_verbosity);

    local_matrix_ = owner_space_->getLocalSpace()->MakeNewSymTMatrix();
  }

  ParSymMatrix::~ParSymMatrix()
  {
    DBG_START_METH("ParSymMatrix::~ParSymMatrix()", dbg_verbosity);
  }

  // assume x is parallel, and y is parallel
  void ParSymMatrix::MultVectorImpl(Number alpha, const Vector &x,
				    Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));

    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

    SmartPtr<const DenseVector> dense_x = par_x->GlobalVector();
    SmartPtr<DenseVector> dense_y = par_y->MakeNewGlobalVector();

    local_matrix_->MultVector(alpha, *dense_x, 0., *dense_y);

    Number *yvalues = dense_y->Values();
    MPI_Allreduce(MPI_IN_PLACE, yvalues, NCols(),
		  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (beta==.0) {
      par_y->ExtractLocalVector(*dense_y);
    }
    else {
      SmartPtr<ParVector> par_v = par_y->MakeNewParVector();
      par_v->ExtractLocalVector(*dense_y);
      par_y->Axpy(beta, *par_v);
    }
  }

  bool ParSymMatrix::HasValidNumbersImpl() const
  {
    int valid = local_matrix_->HasValidNumbers();
    MPI_Allreduce(MPI_IN_PLACE, &valid, 1,
		  MPI_INT, MPI_MIN, MPI_COMM_WORLD);

    return (bool)valid;
  }

  void ParSymMatrix::PrintImpl(const Journalist& jnlst,
                                 EJournalLevel level,
                                 EJournalCategory category,
                                 const std::string& name,
                                 Index indent,
                                 const std::string& prefix) const
  {
    jnlst.PrintfIndented(level, category, indent,
			 "%sParMatrix \"%s\" with %d pieces, nrows %d, ncols:\n",
			 prefix.c_str(), name.c_str(), NumProc(), NRows(), NCols());
    char buffer[256];
    snprintf (buffer, 255, "%s[%2d]", name.c_str(), Rank());
    std::string myname = buffer;
    
    jnlst.StartDistributedOutput();
    local_matrix_->Print(jnlst, level, category, myname, indent+1, prefix);
    jnlst.FinishDistributedOutput();
  }

  ParSymMatrixSpace::ParSymMatrixSpace(Index dim, Index nonZeros, const Index* iRows,
				       const Index* jCols)
    :
    SymMatrixSpace(dim)
  {
    local_space_ = new SymTMatrixSpace(dim, nonZeros, iRows, jCols);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc_);
  }

} // namespace Ipopt
