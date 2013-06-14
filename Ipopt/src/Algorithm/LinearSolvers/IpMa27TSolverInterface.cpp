// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#include "IpoptConfig.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we do not have MA27 in HSL or the linear solver loader, then we want to build the MA27 interface
#if defined(COINHSL_HAS_MA27) || defined(HAVE_LINEARSOLVERLOADER)

#include "IpMa27TSolverInterface.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
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
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

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
    roptions->AddBoundedNumberOption(
      "ma27_pivtol",
      "Pivot tolerance for the linear solver MA27.",
      0.0, true, 1.0, true, 1e-8,
      "A smaller number pivots for sparsity, a larger number pivots for "
      "stability.  This option is only available if Ipopt has been compiled "
      "with MA27.");
    roptions->AddBoundedNumberOption(
      "ma27_pivtolmax",
      "Maximum pivot tolerance for the linear solver MA27.",
      0.0, true, 1.0, true, 1e-4,
      "Ipopt may increase pivtol as high as pivtolmax to get a more accurate "
      "solution to the linear system.  This option is only available if "
      "Ipopt has been compiled with MA27.");
    roptions->AddLowerBoundedNumberOption(
      "ma27_liw_init_factor",
      "Integer workspace memory for MA27.",
      1.0, false, 5.0,
      "The initial integer workspace memory = liw_init_factor * memory "
      "required by unfactored system. Ipopt will increase the workspace "
      "size by meminc_factor if required.  This option is only available if "
      "Ipopt has been compiled with MA27.");
    roptions->AddLowerBoundedNumberOption(
      "ma27_la_init_factor",
      "Real workspace memory for MA27.",
      1.0, false, 5.0,
      "The initial real workspace memory = la_init_factor * memory "
      "required by unfactored system. Ipopt will increase the workspace"
      " size by meminc_factor if required.  This option is only available if "
      " Ipopt has been compiled with MA27.");
    roptions->AddLowerBoundedNumberOption(
      "ma27_meminc_factor",
      "Increment factor for workspace size for MA27.",
      1.0, false, 2.0,
      "If the integer or real workspace is not large enough, "
      "Ipopt will increase its size by this factor.  This option is only "
      "available if Ipopt has been compiled with MA27.");
    roptions->AddStringOption2(
      "ma27_skip_inertia_check",
      "Always pretend inertia is correct.",
      "no",
      "no", "check inertia",
      "yes", "skip inertia check",
      "Setting this option to \"yes\" essentially disables inertia check. "
      "This option makes the algorithm non-robust and easily fail, but it "
      "might give some insight into the necessity of inertia control.");
    roptions->AddStringOption2(
      "ma27_ignore_singularity",
      "Enables MA27's ability to solve a linear system even if the matrix is singular.",
      "no",
      "no", "Don't have MA27 solve singular systems",
      "yes", "Have MA27 solve singular systems",
      "Setting this option to \"yes\" means that Ipopt will call MA27 to "
      "compute solutions for right hand sides, even if MA27 has detected that "
      "the matrix is singular (but is still able to solve the linear system). "
      "In some cases this might be better than using Ipopt's heuristic of "
      "small perturbation of the lower diagonal of the KKT matrix.");

  }

  bool Ma27TSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("ma27_pivtol", pivtol_, prefix);
    if (options.GetNumericValue("ma27_pivtolmax", pivtolmax_, prefix)) {
      ASSERT_EXCEPTION(pivtolmax_>=pivtol_, OPTION_INVALID,
                       "Option \"ma27_pivtolmax\": This value must be between "
                       "ma27_pivtol and 1.");
    }
    else {
      pivtolmax_ = Max(pivtolmax_, pivtol_);
    }

    options.GetNumericValue("ma27_liw_init_factor", liw_init_factor_, prefix);
    options.GetNumericValue("ma27_la_init_factor", la_init_factor_, prefix);
    options.GetNumericValue("ma27_meminc_factor", meminc_factor_, prefix);
    options.GetBoolValue("ma27_skip_inertia_check",
                         skip_inertia_check_, prefix);
    options.GetBoolValue("ma27_ignore_singularity",
                         ignore_singularity_, prefix);
    // The following option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    /* Set the default options for MA27 */
    F77_FUNC(ma27id,MA27ID)(icntl_, cntl_);
#if COIN_IPOPT_VERBOSITY == 0

    icntl_[0] = 0;       // Suppress error messages
    icntl_[1] = 0;       // Suppress diagnostic messages
#endif

    // Reset all private data
    initialized_=false;
    pivtol_changed_ = false;
    refactorize_ = false;

    la_increase_=false;
    liw_increase_=false;

    if (!warm_start_same_structure_) {
      dim_=0;
      nonzeros_=0;
    }
    else {
      ASSERT_EXCEPTION(dim_>0 && nonzeros_>0, INVALID_WARMSTART,
                       "Ma27TSolverInterface called with warm_start_same_structure, but the problem is solved for the first time.");
    }

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
      a_ = NULL;
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

    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    if (!warm_start_same_structure_) {
      dim_ = dim;
      nonzeros_ = nonzeros;

      // Do the symbolic facotrization
      retval = SymbolicFactorization(airn, ajcn);
      if (retval != SYMSOLVER_SUCCESS ) {
        return retval;
      }
    }
    else {
      ASSERT_EXCEPTION(dim_==dim && nonzeros_==nonzeros, INVALID_WARMSTART,
                       "Ma27TSolverInterface called with warm_start_same_structure, but the problem size has changed.");
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus Ma27TSolverInterface::SymbolicFactorization(const Index* airn,
      const Index* ajcn)
  {
    DBG_START_METH("Ma27TSolverInterface::SymbolicFactorization",dbg_verbosity);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // Get memory for the IW workspace
    delete [] iw_;
    iw_ = NULL;

    // Overestimation factor for LIW (20% recommended in MA27 documentation)
    const double LiwFact = 2.0;   // This is 100% overestimation
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "In Ma27TSolverInterface::InitializeStructure: Using overestimation factor LiwFact = %e\n",
                   LiwFact);
    liw_ = (ipfint)(LiwFact*(double(2*nonzeros_+3*dim_+1)));
    iw_ = new ipfint[liw_];

    // Get memory for IKEEP
    delete [] ikeep_;
    ikeep_ = NULL;
    ikeep_ = new ipfint[3*dim_];

    if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                     "\nMatrix structure given to MA27 with dimension %d and %d nonzero entries:\n", dim_, nonzeros_);
      for (Index i=0; i<nonzeros_; i++) {
        Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "A[%5d,%5d]\n",
                       airn[i], ajcn[i]);
      }
    }

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
    const ipfint &iflag = INFO[0];   // Information flag
    const ipfint &ierror = INFO[1];  // Error flag
    const ipfint &nrlnec = INFO[4];  // recommended value for la
    const ipfint &nirnec = INFO[5];  // recommended value for liw

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Return values from MA27AD: IFLAG = %d, IERROR = %d\n",
                   iflag, ierror);

    // Check if error occurred
    if (iflag!=0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "*** Error from MA27AD *** IFLAG = %d IERROR = %d\n", iflag, ierror);
      if (iflag==1) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "The index of a matrix is out of range.\nPlease check your implementation of the Jacobian and Hessian matrices.\n");
      }
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
    }

    // ToDo: try and catch
    // Reserve memory for iw_ for later calls, based on suggested size
    delete [] iw_;
    iw_ = NULL;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Size of integer work space recommended by MA27 is %d\n",
                   nirnec);
    liw_ = (ipfint)(liw_init_factor_ * (double)(nirnec));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Setting integer work space size to %d\n", liw_);
    iw_ = new ipfint[liw_];

    // Reserve memory for a_
    delete [] a_;
    a_ = NULL;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Size of doublespace recommended by MA27 is %d\n",
                   nrlnec);
    la_ = Max(nonzeros_,(ipfint)(la_init_factor_ * (double)(nrlnec)));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Setting double work space size to %d\n", la_);
    a_ = new double[la_];

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

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
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().Start();
    }
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
      iw_ = NULL;
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
    const ipfint &ierror = INFO[1];  // Error flag
    ncmpbr = INFO[11];      // Number of double compressions
    ncmpbi = INFO[12];      // Number of integer compressions
    negevals_ = INFO[14];   // Number of negative eigenvalues

    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "Return values from MA27BD: IFLAG = %d, IERROR = %d\n",
                   iflag, ierror);

    DBG_PRINT((1,"Return from MA27BD iflag = %d and ierror = %d\n",
               iflag, ierror));

    // Check if factorization failed due to insufficient memory space
    // iflag==-3 if LIW too small (recommended value in ierror)
    // iflag==-4 if LA too small (recommended value in ierror)
    if (iflag==-3 || iflag==-4) {
      // Increase size of both LIW and LA
      delete [] iw_;
      iw_ = NULL;
      delete [] a_;
      a_ = NULL;
      ipfint liw_old = liw_;
      ipfint la_old = la_;
      if (iflag==-3) {
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
                     "MA27BD returned iflag=%d and requires more memory.\n Increase liw from %d to %d and la from %d to %d and factorize again.\n",
                     iflag, liw_old, liw_, la_old, la_);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_CALL_AGAIN;
    }

    // Check if the system is singular, and if some other error occurred
    if (iflag==-5 || (!ignore_singularity_ && iflag==3)) {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_SINGULAR;
    }
    else if (iflag==3) {
      Index missing_rank = dim_-INFO[1];
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned iflag=%d and detected rank deficiency of degree %d.\n",
                     iflag, missing_rank);
      // We correct the number of negative eigenvalues here to include
      // the zero eigenvalues, since otherwise we indicate the wrong
      // inertia.
      negevals_ += missing_rank;
    }
    else if (iflag != 0) {
      // There is some error
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
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

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of doubles for MA27 to hold factorization (INFO(9)) = %d\n",
                   INFO[8]);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Number of integers for MA27 to hold factorization (INFO(10)) = %d\n",
                   INFO[9]);

    // Check whether the number of negative eigenvalues matches the requested
    // count
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemFactorization().End();
    }
    if (!skip_inertia_check_ && check_NegEVals && (numberOfNegEVals!=negevals_)) {
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
    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().Start();
    }

    ipfint N=dim_;
    double* W = new double[maxfrt_];
    ipfint* IW1 = new ipfint[nsteps_];

    // For each right hand side, call MA27CD
    for (Index irhs=0; irhs<nrhs; irhs++) {
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

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemBackSolve().End();
    }
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
                   "Indreasing pivot tolerance for MA27 from %7.2e ",
                   pivtol_);
    pivtol_ = Min(pivtolmax_, pow(pivtol_,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   pivtol_);
    return true;
  }

} // namespace Ipopt

#endif /* COINHSL_HAS_MA27 or HAVE_LINEARSOLVERLOADER */
