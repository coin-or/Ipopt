// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Sanjeeb Dash, Andreas Waechter     IBM    2009-06-17

#include "IpParExpansionMatrix.hpp"
#include "IpParVector.hpp"

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

    local_matrix_ = owner_space_->getLocalSpace()->MakeNewExpansionMatrix();
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

  }

  void ParExpansionMatrix::TransMultVectorImpl(Number alpha, const Vector &x,
      Number beta, Vector &y) const
  {
    const ParVector* par_x = static_cast<const ParVector*>(&x);
    DBG_ASSERT(dynamic_cast<const ParVector*>(&x));
    ParVector* par_y = static_cast<ParVector*>(&y);
    DBG_ASSERT(dynamic_cast<ParVector*>(&y));

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

  }

  void ParExpansionMatrix::ComputeRowAMaxImpl(Vector& rows_norms, bool init) const
  {
    ParVector* par_vec = static_cast<ParVector*>(&rows_norms);
    DBG_ASSERT(dynamic_cast<ParVector*>(&rows_norms));
    Number* vec_vals=par_vec->LocalVector()->Values();

  }

  void ParExpansionMatrix::ComputeColAMaxImpl(Vector& cols_norms, bool init) const
  {
  }

  void ParExpansionMatrix::PrintImpl(const Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category,
                                  const std::string& name,
                                  Index indent,
                                  const std::string& prefix) const
  {
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sParExpansionMatrix \"%s\" with %d rows and %d columns:\n",
                         prefix.c_str(), name.c_str(), NRows(), NCols());

    for (Index i=0; i<NCols(); i++) {
      jnlst.PrintfIndented(level, category, indent,
                           "%s%s[%5d,%5d]=%23.16e  (%d)\n",
                           prefix.c_str(), name.c_str(), /*exp_pos[i]+*/1,
                           i+1, 1., i);
    }
  }

  ParExpansionMatrixSpace::ParExpansionMatrixSpace(SmartPtr<const ParVectorSpace> LargeVectorSpace, 
						   SmartPtr<const ParVectorSpace> SmallVectorSpace, 
						   const Index *ExpPos,
			  const int offset)
    :
    MatrixSpace(LargeVectorSpace->LocalSize(), SmallVectorSpace->LocalSize())
  {
    local_space_ = new ExpansionMatrixSpace(LargeVectorSpace->LocalSize(),
					    SmallVectorSpace->LocalSize(),
					    ExpPos, offset);

  }

} // namespace Ipopt
