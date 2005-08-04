// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMa27TSolverInterface.hpp"

#ifdef OLD_C_HEADERS
# include <math.h>
#else
# include <cmath>
#endif

/** Prototypes for MA27's Fortran subroutines */
extern "C"
{
  void F77_FUNC(ma27id,MA27ID)(ipfint* ICNTL, double* CNTL);
  void F77_FUNC(ma27ad,MA27AD)(ipfint *N, ipfint *NZ, const ipfint *IRN, const ipfint* ICN,
                               ipfint *IW, ipfint* LIW, ipfint* IKEEP, ipfint *IW1,
                               ipfint* NSTEPS, ipfint* IFLAG, ipfint* ICNTL,
                               double* CNTL, ipfint *INFO, double* OPS);
  void F77_FUNC(ma27bd,MA27BD)(ipfint *N, ipfint *NZ, const ipfint *IRN, const ipfint* ICN,
                               double* A, ipfint* LA, ipfint* IW, ipfint* LIW,
                               ipfint* IKEEP, ipfint* NSTEPS, ipfint* MAXFRT,
                               ipfint* IW1, ipfint* ICNTL, double* CNTL,
                               ipfint* INFO);
  void F77_FUNC(ma27cd,MA27CD)(ipfint *N, double* A, ipfint* LA, ipfint* IW,
                               ipfint* LIW, double* W, ipfint* MAXFRT,
                               double* RHS, ipfint* IW1, ipfint* NSTEPS,
                               ipfint* ICNTL, double* CNTL);
}

namespace Ipopt
{

  DBG_SET_VERBOSITY(0);

  DefineIpoptType(Ma27TSolverInterface);

  Ma27TSolverInterface::Ma27TSolverInterface()
      :
      dim_(0),
      nonzeros_(0),
      initialized_(false),
      pivtol_changed_(false),
      refactorize_(false),

      liw_(0),
      iw_(NULL),
      ikeep_(NULL),
      la_(0),
      a_(NULL),

      la_increase_(false),
      liw_increase_(false)
  {
    DBG_START_METH("Ma27TSolverInterface::Ma27TSolverInterface()",dbg_verbosity);
  }

  Ma27TSolverInterface::~Ma27TSolverInterface()
  {
    DBG_START_METH("Ma27TSolverInterface::~Ma27TSolverInterface()",
                   dbg_verbosity);
    delete [] iw_;
    delete [] ikeep_;
    delete [] a_;
  }

  void Ma27TSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption("pivtol", "pivot tolerance for the linear solver. smaller number - pivot for sparsity, larger number - pivot for stability",
                                     0.0, true, 1.0, true, 1e-8);
    roptions->AddBoundedNumberOption("pivtolmax", "maximum pivot tolerance. IPOPT may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system",
                                     0.0, true, 1.0, true, 1e-4);
    roptions->AddLowerBoundedNumberOption("liw_init_factor", "integer workspace memory = liw_init_factor * memory required by unfactored system",
                                          1.0, false, 5.0);
    roptions->AddLowerBoundedNumberOption("la_init_factor", "real workspace memory = la_init_factor * memory required by unfactored system",
                                          1.0, false, 5.0);
    roptions->AddLowerBoundedNumberOption("meminc_factor", "if workspace is not large enough, IPOPT will increase the size by this factor",
                                          1.0, false, 10.0);
  }

  bool Ma27TSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("pivtol", pivtol_, prefix);
    if(options.GetNumericValue("pivtolmax", pivtolmax_, prefix)) {
      ASSERT_EXCEPTION(pivtolmax_>=pivtol_, OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"pivtolmax\": This value must be between pivtol and 1.");
    }
    else {
      pivtolmax_ = Max(pivtolmax_, pivtol_);
    }

    options.GetNumericValue("liw_init_factor", liw_init_factor_, prefix);
    options.GetNumericValue("la_init_factor", la_init_factor_, prefix);
    options.GetNumericValue("meminc_factor", meminc_factor_, prefix);

    /* Set the default options for MA27 */
    F77_FUNC(ma27id,MA27ID)(icntl_, cntl_);
#ifndef IP_DEBUG

    icntl_[0] = 0;       // Suppress error messages
    icntl_[1] = 0;       // Suppress diagnostic messages
#endif

    // Reset all private data
    dim_=0;
    nonzeros_=0;
    initialized_=false;
    pivtol_changed_ = false;
    refactorize_ = false;

    la_increase_=false;
    liw_increase_=false;

    return true;
  }

  ESymSolverStatus Ma27TSolverInterface::MultiSolve(bool new_matrix,
      const Index* airn,
      const Index* ajcn,
      Index nrhs,
      double* rhs_vals,
      bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("Ma27TSolverInterface::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    DBG_ASSERT(initialized_);
    DBG_ASSERT(la_!=0);

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
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(airn, ajcn, check_NegEVals, numberOfNegEVals);
      if (retval!=SYMSOLVER_SUCCESS) {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
        return retval;  // Matrix singular or error occurred
      }
      refactorize_ = false;
    }

    // do the backsolve
    return Backsolve(nrhs, rhs_vals);
  }

  double* Ma27TSolverInterface::GetValuesArrayPtr()
  {
    DBG_START_METH("Ma27TSolverInterface::GetValuesArrayPtr",dbg_verbosity);
    DBG_ASSERT(initialized_);

    // If the size of a is to be increase for the next factorization
    // anyway, delete the current large array and just return enough
    // to store the values

    if (la_increase_) {
      delete [] a_;
      a_ = new double [nonzeros_];
    }

    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus Ma27TSolverInterface::InitializeStructure(Index dim, Index nonzeros,
      const Index* airn,
      const Index* ajcn)
  {
    DBG_START_METH("Ma27TSolverInterface::InitializeStructure",dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Do the symbolic facotrization
    ESymSolverStatus retval = SymbolicFactorization(airn, ajcn);
    if (retval != SYMSOLVER_SUCCESS ) {
      return retval;
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus Ma27TSolverInterface::SymbolicFactorization(const Index* airn,
      const Index* ajcn)
  {
    DBG_START_METH("Ma27TSolverInterface::SymbolicFactorization",dbg_verbosity);

    // Get memory for the IW workspace
    delete [] iw_;

    // Overstimation factor for LIW (20% recommended in MA27 documentation)
    const double LiwFact = 2.0;   // This is 100% overestimation
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "In Ma27TSolverInterface::InitializeStructure: Using overestimation factor LiwFact = %e\n",
                   LiwFact);
    liw_ = (ipfint)(LiwFact*(double(2*nonzeros_+3*dim_+1)));
    iw_ = new ipfint[liw_];

    // Get memory for IKEEP
    delete [] ikeep_;
    ikeep_ = new ipfint[3*dim_];

    // Call MA27AD (cast to ipfint for Index types)
    ipfint N = dim_;
    ipfint NZ = nonzeros_;
    ipfint IFLAG = 0;
    double OPS;
    ipfint INFO[20];
    ipfint* IW1 = new ipfint[2*dim_];  // Get memory for IW1 (only local)
    F77_FUNC(ma27ad,MA27AD)(&N, &NZ, airn, ajcn, iw_, &liw_, ikeep_,
                            IW1, &nsteps_, &IFLAG, icntl_, cntl_,
                            INFO, &OPS);
    delete [] IW1;  // No longer required

    // Receive several information
    ipfint iflag = INFO[0];   // Information flag
    ipfint ierror = INFO[1];  // Error flag
    ipfint nrlnec = INFO[4];  // recommended value for la
    ipfint nirnec = INFO[5];  // recommended value for liw

    // Check if error occurred
    if (iflag!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "*** Error from MA27AD *** IFLAG = %d IERROR = %d\n", iflag, ierror);
      return SYMSOLVER_FATAL_ERROR;
    }

    // ToDo: try and catch
    // Reserve memory for iw_ for later calls, based on suggested size
    delete [] iw_;
    liw_ = (ipfint)(liw_init_factor_ * (double)(nirnec));
    iw_ = new ipfint[liw_];

    // Reserve memory for a_
    delete [] a_;
    la_ = Max(nonzeros_,(ipfint)(la_init_factor_ * (double)(nrlnec)));
    a_ = new double[la_];

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  Ma27TSolverInterface::Factorization(const Index* airn,
                                      const Index* ajcn,
                                      bool check_NegEVals,
                                      Index numberOfNegEVals)
  {
    DBG_START_METH("Ma27TSolverInterface::Factorization",dbg_verbosity);
    // Check if la should be increased
    if (la_increase_) {
      double* a_old = a_;
      ipfint la_old = la_;
      la_ = (ipfint)(meminc_factor_ * (double)(la_));
      a_ = new double[la_];
      for (Index i=0; i<nonzeros_; i++) {
        a_[i] = a_old[i];
      }
      delete [] a_old;
      la_increase_ = false;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: Increasing la from %d to %d\n",
                     la_old, la_);
    }

    // Check if liw should be increased
    if (liw_increase_) {
      delete [] iw_;
      ipfint liw_old = liw_;
      liw_ = (ipfint)(meminc_factor_ * (double)(liw_));
      iw_ = new ipfint[liw_];
      liw_increase_ = false;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: Increasing liw from %d to %d\n",
                     liw_old, liw_);
    }

    ipfint iflag;  // Information flag
    ipfint ncmpbr; // Number of double precision compressions
    ipfint ncmpbi; // Number of integer compressions

    // Call MA27BD; possibly repeatedly if workspaces are too small
    ipfint N=dim_;
    ipfint NZ=nonzeros_;
    ipfint* IW1 = new ipfint[2*dim_];
    ipfint INFO[20];
    cntl_[0] = pivtol_;  // Set pivot tolerance

    F77_FUNC(ma27bd,MA27BD)(&N, &NZ, airn, ajcn, a_,
                            &la_, iw_, &liw_, ikeep_, &nsteps_,
                            &maxfrt_, IW1, icntl_, cntl_, INFO);
    delete [] IW1;

    // Receive information about the factorization
    iflag = INFO[0];        // Information flag
    ipfint ierror = INFO[1];  // Error flag
    ncmpbr = INFO[11];      // Number of double compressions
    ncmpbi = INFO[12];      // Number of integer compressions
    negevals_ = INFO[14];   // Number of negative eigenvalues

    DBG_PRINT((1,"Return from MA27 iflag = %d and ierror = %d\n",
               iflag, ierror));

    // Check if factorization failed due to insufficient memory space
    // iflag==-3 if LIW too small (recommended value in ierror)
    // iflag==-4 if LA too small (recommended value in ierror)
    if (iflag==-3 || iflag==-4) {
      // Increase size of both LIW and LA
      delete [] iw_;
      delete [] a_;
      ipfint liw_old = liw_;
      ipfint la_old = la_;
      if(iflag==-3) {
        liw_ = (ipfint)(meminc_factor_ * (double)(ierror));
        la_ = (ipfint)(meminc_factor_ * (double)(la_));
      }
      else {
        liw_ = (ipfint)(meminc_factor_ * (double)(liw_));
        la_ = (ipfint)(meminc_factor_ * (double)(ierror));
      }
      iw_ = new ipfint[liw_];
      a_ = new double[la_];
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned iflag=%d.\n Increase liw from %d to %d and la from %d to %d and factorize again.\n",
                     iflag, liw_old, liw_, la_old, la_);
      return SYMSOLVER_CALL_AGAIN;
    }

    // Check if the system is singular, and if some other error occurred
    if (iflag==-5 || iflag==3) {
      return SYMSOLVER_SINGULAR;
    }
    else if (iflag != 0) {
      // There is some error
      return SYMSOLVER_FATAL_ERROR;
    }

    // Check if it might be more efficient to use more memory next time
    // (if there were too many compressions for this factorization)
    if (ncmpbr>=10) {
      la_increase_ = true;
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned ncmpbr=%d. Increase la before the next factorization.\n",
                     ncmpbr);
    }
    if (ncmpbi>=10) {
      liw_increase_ = true;
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned ncmpbi=%d. Increase liw before the next factorization.\n",
                     ncmpbr);
    }

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: negevals_ = %d, but numberOfNegEVals = %d\n",
                     negevals_, numberOfNegEVals);
      return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus Ma27TSolverInterface::Backsolve(Index nrhs,
      double *rhs_vals)
  {
    DBG_START_METH("Ma27TSolverInterface::Backsolve",dbg_verbosity);

    ipfint N=dim_;
    double* W = new double[maxfrt_];
    ipfint* IW1 = new ipfint[nsteps_];

    // For each right hand side, call MA27CD
    for(Index irhs=0; irhs<nrhs; irhs++) {
      if (DBG_VERBOSITY()>=2) {
        for (Index i=0; i<dim_; i++) {
          DBG_PRINT((2, "rhs[%5d] = %23.15e\n", i, rhs_vals[irhs*dim_+i]));
        }
      }

      F77_FUNC(ma27cd,MA27CD)(&N, a_, &la_, iw_, &liw_, W, &maxfrt_,
                              &rhs_vals[irhs*dim_], IW1, &nsteps_,
                              icntl_, cntl_);

      if (DBG_VERBOSITY()>=2) {
        for (Index i=0; i<dim_; i++) {
          DBG_PRINT((2, "sol[%5d] = %23.15e\n", i, rhs_vals[irhs*dim_+i]));
        }
      }
    }
    delete [] W;
    delete [] IW1;

    return SYMSOLVER_SUCCESS;
  }

  Index Ma27TSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("Ma27TSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(ProvidesInertia());
    DBG_ASSERT(initialized_);
    return negevals_;
  }

  bool Ma27TSolverInterface::IncreaseQuality()
  {
    DBG_START_METH("Ma27TSolverInterface::IncreaseQuality",dbg_verbosity);
    if (pivtol_ == pivtolmax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Indreasing pivot tolerance from %7.2e ",
                   pivtol_);
    pivtol_ = Min(pivtolmax_, pow(pivtol_,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   pivtol_);
    return true;
  }

} // namespace Ipopt
