// Copyright (C) 2006, 2009 Damien Hocking, KBC Advanced Technologies and others
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter, Sanjeeb Dash     IBM    2009-08-11
//                  (based on IpParMumpsSolverInterface.hpp rev 1543)

#include "IpParGenMumpsSolverInterface.hpp"

#include "IpMpi.hpp"

#include "dmumps_c.h"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#ifdef HAVE_CSTDLIB
# include <cstdlib>
#else
# ifdef HAVE_STDLIB_H
#  include <stdlib.h>
# else
#  error "don't have header file for stdlib"
# endif
#endif

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  ParGenMumpsSolverInterface::ParGenMumpsSolverInterface()
  {
    DBG_START_METH("ParGenMumpsSolverInterface::ParGenMumpsSolverInterface()",
                   dbg_verbosity);
    //initialize mumps
    DMUMPS_STRUC_C* mumps_ = new DMUMPS_STRUC_C;
    mumps_->n = 0;
    mumps_->nz_loc = 0;
    mumps_->a_loc = NULL;
    mumps_->jcn_loc = NULL;
    mumps_->irn_loc = NULL;
    mumps_->job = -1;//initialize mumps
    mumps_->par = 1;//working host for sequential version
    mumps_->sym = 0;//unsymmetric matrix!!!
    mumps_->comm_fortran = -987654; // This prompts MUMPS to use MPI_COMM_WORLD
    dmumps_c(mumps_);
    mumps_ptr_ = (void*)mumps_;
  }


  ParGenMumpsSolverInterface::~ParGenMumpsSolverInterface()
  {
    DBG_START_METH("ParGenMumpsSolverInterface::~ParGenMumpsSolverInterface()",
                   dbg_verbosity);

    DMUMPS_STRUC_C* mumps_ = (DMUMPS_STRUC_C*)mumps_ptr_;
    mumps_->job = -2; //terminate mumps
    dmumps_c(mumps_);
    delete [] mumps_->a_loc;
    delete mumps_;
  }

  void ParGenMumpsSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool ParGenMumpsSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // all options defined in MumpsSolverInterface
    options.GetNumericValue("mumps_pivtol", pivtol_, prefix);
    if (options.GetNumericValue("mumps_pivtolmax", pivtolmax_, prefix)) {
      ASSERT_EXCEPTION(pivtolmax_>=pivtol_, OPTION_INVALID,
                       "Option \"mumps_pivtolmax\": This value must be between "
                       "mumps_pivtol and 1.");
    }
    else {
      pivtolmax_ = Max(pivtolmax_, pivtol_);
    }

    options.GetIntegerValue("mumps_mem_percent",
                            mem_percent_, prefix);

    // The following option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    options.GetIntegerValue("mumps_permuting_scaling",
                            mumps_permuting_scaling_, prefix);
    options.GetIntegerValue("mumps_pivot_order", mumps_pivot_order_, prefix);
    options.GetIntegerValue("mumps_scaling", mumps_scaling_, prefix);
    options.GetNumericValue("mumps_dep_tol", mumps_dep_tol_, prefix);

    // Reset all private data
    initialized_ = false;
    pivtol_changed_ = false;
    refactorize_ = false;
    have_symbolic_factorization_ = false;

    DMUMPS_STRUC_C* mumps_ = (DMUMPS_STRUC_C*)mumps_ptr_;
    if (!warm_start_same_structure_) {
      mumps_->n = 0;
      mumps_->nz_loc = 0;
    }
    else {
      ASSERT_EXCEPTION(mumps_->n>0 && mumps_->nz_loc>0, INVALID_WARMSTART,
                       "ParGenMumpsSolverInterface called with warm_start_same_structure, but the problem is solved for the first time.");
    }

    return true;
  }

  ESymSolverStatus ParGenMumpsSolverInterface::MultiSolve(bool new_matrix,
      const Index* ia,
      const Index* ja,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("ParGenMumpsSolverInterface::MultiSolve", dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);
    DBG_ASSERT(((DMUMPS_STRUC_C*)mumps_ptr_)->irn_loc == ia);
    DBG_ASSERT(((DMUMPS_STRUC_C*)mumps_ptr_)->jcn_loc == ja);

    if (pivtol_changed_) {
      DBG_PRINT((1,"Pivot tolerance has changed.\n"));
      pivtol_changed_ = false;
      // If the pivot tolerance has been changed but the matrix is not
      // new, we have to request the values for the matrix again to do
      // the factorization again.
      if (!new_matrix) {
        DBG_PRINT((1,"Ask caller to call again.\n"));
        refactorize_ = true;
        return SYMSOLVER_CALL_AGAIN;
      }
    }

    // check if a factorization has to be done
    DBG_PRINT((1, "new_matrix = %d\n", new_matrix));
    if (new_matrix || refactorize_) {
      ESymSolverStatus retval;
      // Do the symbolic facotrization if it hasn't been done yet
      if (!have_symbolic_factorization_) {
        retval = SymbolicFactorization();
        if (retval != SYMSOLVER_SUCCESS ) {
          return retval;
        }
        have_symbolic_factorization_ = true;
      }
      // perform the factorization
      retval = Factorization(check_NegEVals, numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS)  {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
        return retval;  // Matrix singular or error occurred
      }
      refactorize_ = false;
    }
    // do the solve
    return Solve(nrhs, rhs_vals);
  }


  double* ParGenMumpsSolverInterface::GetValuesArrayPtr()
  {
    DMUMPS_STRUC_C* mumps_ = (DMUMPS_STRUC_C*)mumps_ptr_;
    DBG_START_METH("ParGenMumpsSolverInterface::GetValuesArrayPtr",dbg_verbosity)
    DBG_ASSERT(initialized_);
    return mumps_->a_loc;
  }

  /* Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus ParGenMumpsSolverInterface::InitializeStructure(Index dim,
      Index nonzeros,
      const Index* ia,
      const Index* ja)
  {
    DMUMPS_STRUC_C* mumps_ = (DMUMPS_STRUC_C*)mumps_ptr_;
    DBG_START_METH("ParGenMumpsSolverInterface::InitializeStructure", dbg_verbosity);

    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    if (!warm_start_same_structure_) {
      mumps_->n = dim;
      mumps_->nz_loc = nonzeros;
      delete [] mumps_->a_loc;
      mumps_->a_loc = NULL;

      mumps_->a_loc = new double[nonzeros];
      mumps_->irn_loc = const_cast<int*>(ia);
      mumps_->jcn_loc = const_cast<int*>(ja);

      // make sure we do the symbolic factorization before a real
      // factorization
      have_symbolic_factorization_ = false;
    }
    else {
      ASSERT_EXCEPTION(mumps_->n==dim && mumps_->nz_loc==nonzeros,
                       INVALID_WARMSTART,"ParGenMumpsSolverInterface called with warm_start_same_structure, but the problem size has changed.");
    }

    initialized_ = true;
    return retval;
  }


  ESymSolverStatus ParGenMumpsSolverInterface::SymbolicFactorization()
  {
    DBG_START_METH("ParGenMumpsSolverInterface::SymbolicFactorization",
                   dbg_verbosity);
    DMUMPS_STRUC_C* mumps_data = (DMUMPS_STRUC_C*)mumps_ptr_;

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    mumps_data->job = 1;//symbolic ordering pass

    mumps_data->icntl[1] = 0;
    mumps_data->icntl[2] = 0;//QUIETLY!
    mumps_data->icntl[3] = 0;

    //mumps_data->icntl[0] = 6;
    //mumps_data->icntl[1] = 6;
    //mumps_data->icntl[2] = 6;//QUIETLY!
    //mumps_data->icntl[3] = 4;

    mumps_data->icntl[5] = mumps_permuting_scaling_;
    mumps_data->icntl[6] = mumps_pivot_order_;
    mumps_data->icntl[7] = mumps_scaling_;
    // TODO: WHAT IS GOOD FOR UNSYMMETIC SYSTEMS?
    mumps_data->icntl[9] = 0;//no iterative refinement iterations


    mumps_data->icntl[12] = 0;//we don't need the inertia and can use Scalapack
    mumps_data->icntl[13] = mem_percent_; //% memory to allocate over expected

    mumps_data->icntl[17] = 3; //matrix is provided distributedly

    mumps_data->cntl[0] = pivtol_;  // Set pivot tolerance

    // not implemented in parallel
    // dump_matrix(mumps_data);

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling MUMPS-1 for symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    dmumps_c(mumps_data);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with MUMPS-1 for symbolic factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    int error = mumps_data->infog[0];
    const int& mumps_permuting_scaling_used = mumps_data->infog[22];
    const int& mumps_pivot_order_used = mumps_data->infog[6];
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "MUMPS used permuting_scaling %d and pivot_order %d.\n",
                   mumps_permuting_scaling_used, mumps_pivot_order_used);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "           scaling will be %d.\n",
                   mumps_data->icntl[7]);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

    //return appropriat value
    if (error == -6) {//system is singular
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "MUMPS returned INFOG(1) = %d matrix is singular.\n",error);
      return SYMSOLVER_SINGULAR;
    }
    if (error < 0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error=%d returned from MUMPS in Factorization.\n",
                     error);
      return SYMSOLVER_FATAL_ERROR;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus ParGenMumpsSolverInterface::Factorization(
    bool check_NegEVals, Index numberOfNegEVals)
  {
    DBG_START_METH("ParGenMumpsSolverInterface::Factorization", dbg_verbosity);
    DMUMPS_STRUC_C* mumps_data = (DMUMPS_STRUC_C*)mumps_ptr_;
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().Start();
    }

    mumps_data->job = 2;//numerical factorization

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Calling MUMPS-2 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    dmumps_c(mumps_data);
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Done with MUMPS-2 for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
    int error = mumps_data->infog[0];

    //Check for errors
    if (error == -8 || error == -9) {//not enough memory
      const Index trycount_max = 20;
      for (int trycount=0; trycount<trycount_max; trycount++) {
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                       "MUMPS returned INFO(1) = %d and requires more memory, reallocating.  Attempt %d\n",
                       error,trycount+1);
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                       "  Increasing icntl[13] from %d to ", mumps_data->icntl[13]);
        double mem_percent = mumps_data->icntl[13];
        mumps_data->icntl[13] = (Index)(2.0 * mem_percent);
        Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA, "%d.\n", mumps_data->icntl[13]);

        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "Calling MUMPS-2 (repeated) for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
        dmumps_c(mumps_data);
        Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                       "Done with MUMPS-2 (repeated) for numerical factorization at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
        error = mumps_data->infog[0];
        if (error != -8 && error != -9)
          break;
      }
      if (error == -8 || error == -9) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "MUMPS was not able to obtain enough memory.\n");
        if (HaveIpData()) {
          IpData().TimingStats().LinearSystemFactorization().End();
        }
        return SYMSOLVER_FATAL_ERROR;
      }
    }
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().End();
    }

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of doubles for MUMPS to hold factorization (INFOG(9)) = %d\n",
                   mumps_data->infog[8]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of integers for MUMPS to hold factorization (INFOG(10)) = %d\n",
                   mumps_data->infog[9]);

    if (error == -10) {//system is singular
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "MUMPS returned INFOG(1) = %d matrix is singular.\n",error);
      return SYMSOLVER_SINGULAR;
    }

    negevals_ = mumps_data->infog[11];

    if (error == -13) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "MUMPS returned INFOG(1) =%d - out or memory.\nIn some cases it helps to decrease the value of the option \"mumps_mem_percent\".\n",
                     error);
      return SYMSOLVER_FATAL_ERROR;
    }
    if (error < 0) {//some other error
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "MUMPS returned INFOG(1) =%d MUMPS failure.\n",
                     error);
      return SYMSOLVER_FATAL_ERROR;
    }

    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In ParGenMumpsSolverInterface::Factorization: negevals_ = %d, but numberOfNegEVals = %d\n",
                     negevals_, numberOfNegEVals);
      return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus ParGenMumpsSolverInterface::Solve(Index nrhs, double *rhs_vals)
  {
    DBG_START_METH("ParGenMumpsSolverInterface::Solve", dbg_verbosity);
    DMUMPS_STRUC_C* mumps_data = (DMUMPS_STRUC_C*)mumps_ptr_;
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }
    for (Index i = 0; i < nrhs; i++) {
      Index offset = i * mumps_data->n;
      mumps_data->rhs = &(rhs_vals[offset]);
      mumps_data->job = 3;//solve
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Calling MUMPS-3 for solve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
      dmumps_c(mumps_data);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Done with MUMPS-3 for solve at cpu time %10.3f (wall %10.3f).\n", CpuTime(), WallclockTime());
      int error = mumps_data->infog[0];
      if (error < 0) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error=%d returned from MUMPS in Solve.\n",
                       error);
        retval = SYMSOLVER_FATAL_ERROR;
      }
      else {
        // distribute data to all processors
        MPI_Bcast(&rhs_vals[offset], mumps_data->n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      }
    }
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().End();
    }
    return retval;
  }

  bool ParGenMumpsSolverInterface::IncreaseQuality()
  {
    DBG_START_METH("MumpsTSolverInterface::IncreaseQuality",dbg_verbosity);
    if (pivtol_ == pivtolmax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Increasing pivot tolerance for MUMPS from %7.2e ",
                   pivtol_);

    //this is a more aggresive update then MA27
    //this should be tuned
    pivtol_ = Min(pivtolmax_, pow(pivtol_,0.5));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   pivtol_);
    return true;
  }

}//end Ipopt namespace



