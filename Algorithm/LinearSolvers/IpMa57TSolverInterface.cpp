// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Michael Hagemann               Univ of Basel 2005-10-28
//               original version (based on MA27TSolverInterface.cpp)

#include "IpMa57TSolverInterface.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <iostream>


/** Prototypes for MA57's Fortran subroutines */
extern "C"
{
  /*
   *  MA57ID -- Initialize solver.
   */
  extern void  F77_FUNC (ma57id, MA57ID) (
      double    *cntl,
      ipfint    *icntl);

  /*
   *  MA57AD -- Symbolic Factorization.
   */
  extern void  F77_FUNC (ma57ad, MA57AD) (
      ipfint    *n,     /* Order of matrix. */
      ipfint    *ne,            /* Number of entries. */

      const ipfint    *irn,       /* Matrix nonzero row structure */
      const ipfint    *jcn,       /* Matrix nonzero column structure */

      ipfint    *lkeep,     /* Workspace for the pivot order of lenght 3*n */
      ipfint    *keep,      /* Workspace for the pivot order of lenght 3*n */
      /* Automatically iflag = 0; ikeep pivot order iflag = 1 */
      ipfint    *iwork,     /* Integer work space. */
      ipfint    *icntl,     /* Integer Control parameter of length 30*/
      ipfint    *info,      /* Statistical Information; Integer array of length 20 */
      double    *rinfo);    /* Double Control parameter of length 5 */

  /*
   * MA57BD -- Numerical Factorization.
   */
  extern void  F77_FUNC (ma57bd, MA57BD) (
      ipfint    *n,     /* Order of matrix. */
      ipfint    *ne,            /* Number of entries. */

      double    *a,     /* Numerical values. */
      double    *fact,      /* Entries of factors. */
      ipfint    *lfact,     /* Length of array `fact'. */
      ipfint    *ifact,     /* Indexing info for factors. */
      ipfint    *lifact,    /* Length of array `ifact'. */

      ipfint    *lkeep,     /* Length of array `keep'. */
      ipfint    *keep,      /* Integer array. */

      ipfint    *iwork,     /* Workspace of length `n'. */

      ipfint    *icntl,     /* Integer Control parameter of length 20. */
      double    *cntl,      /* Double Control parameter of length 5. */
      ipfint    *info,      /* Statistical Information; Integer array of length 40. */
      double    *rinfo);    /* Statistical Information; Real array of length 20. */

  /*
   * MA57CD -- Solution.
   */
  extern void  F77_FUNC (ma57cd, MA57CD) (
      ipfint    *job,       /* Solution job.  Solve for... */
      /* JOB <= 1:  A */
      /* JOB == 2:  PLP^t */
      /* JOB == 3:  PDP^t */
      /* JOB >= 4:  PL^t P^t */

      ipfint    *n,         /* Order of matrix. */

      double    *fact,      /* Entries of factors. */
      ipfint    *lfact,     /* Length of array `fact'. */
      ipfint    *ifact,     /* Indexing info for factors. */
      ipfint    *lifact,    /* Length of array `ifact'. */

      ipfint    *nrhs,      /* Number of right hand sides. */
      double    *rhs,       /* Numerical Values. */
      ipfint    *lrhs,      /* Leading dimensions of `rhs'. */

      double    *work,      /* Real workspace. */
      ipfint    *lwork,     /* Length of `work', >= N*NRHS. */
      ipfint    *iwork,     /* Integer array of length `n'. */

      ipfint    *icntl,     /* Integer Control parameter array of length 20. */
      ipfint    *info);     /* Statistical Information; Integer array of length 40. */

  /*
   * MC57ED -- Copy arrays.
   */
  extern void  F77_FUNC (ma57ed, MA57ED) (
      ipfint    *n,
      ipfint    *ic,        /* 0: copy real array.  >=1:  copy integer array. */
      ipfint    *keep,

      double    *fact,
      ipfint    *lfact,
      double    *newfac,
      ipfint    *lnew,

      ipfint    *ifact,
      ipfint    *lifact,
      ipfint    *newifc,
      ipfint    *linew,

      ipfint    *info);
}

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  std::string ma57_err_msg[] = {
                                 "Operation successful.\n",

                                 "Value of N is out of range on a call to MA57A/AD, MA57B/BD, MA57C/CD, or\n"
                                 "MA57D/DD. Value given is held in INFO(2).\n",

                                 "Value of NE is out of range on a call to MA57A/AD, MA57B/BD, or\n"
                                 "MA57D/DD. Value given is held in INFO(2).\n",

                                 "Failure due to insufficient REAL space on a call to MA57B/BD. INFO(17)\n"
                                 "is set to a value that may suffice. INFO(2) is set to value of\n"
                                 "LFACT. The user can allocate a larger array and copy the contents of\n"
                                 "FACT into it using MA57E/ED, and recall MA57B/BD.\n",

                                 "Failure due to insufficient INTEGER space on a call to\n"
                                 "MA57B/BD. INFO(18) is set to a value that may suffice. INFO(2) is set to\n"
                                 "value of LIFACT. The user can allocate a larger array and copy the\n"
                                 "contents of IFACT into it using MA57E/ED, and recall MA57B/BD.\n",

                                 "A pivot with magnitude less than or equal to CNTL(2) was found at pivot\n"
                                 "step INFO(2) when calling MA57B/BD with ICNTL(7) = 2 or 3, or the\n"
                                 "correction obtained when using matrix modification does not give a pivot\n"
                                 "greater than CNTL(2) when ICNTL(7) = 4.\n",

                                 "A change of sign of pivots has been detected when ICNTL(7) = 2. INFO(2)\n"
                                 "is set to the pivot step at which the change was detected on a call to\n"
                                 "MA57B/BD.\n",

                                 "Either LNEW < LFACT or LINEW < LIFACT on a call to MA57E/ED. INFO(2) is\n"
                                 "set to LNEW or LINEW as appropriate.\n",

                                 "Iterative refinement fails to converge in specified number of iterations\n"
                                 "on a call to MA57D/DD.\n",

                                 "Error in permutation array when ICNTL(6)=1 on a call to\n"
                                 "MA57A/AD. INFO(2) holds first component at which error was detected.\n",

                                 "Value of ICNTL(7) out of range on a call to MA57B/BD. Value given held\n"
                                 "in INFO(2).\n",

                                 "LRHS < N on a call to MA57C/CD. INFO(2) holds value of LRHS.\n",

                                 "Invalid value for JOB on a call to MA57D/DD. Value given held in\n"
                                 "INFO(2).\n",

                                 "Invalid value of ICNTL(9) on a call to MA57D/DD. Value given held in\n"
                                 "INFO(2).\n",

                                 "Failure of MC71A/AD on a call to MA57D/DD with ICNTL(10)> 0.\n",

                                 "LKEEP less than 5*N+NE+MAX(N,NE) +42 on a call to MA57A/AD or\n"
                                 "MA57B/BD. INFO(2) holds value of LKEEP.\n",

                                 "NRHS less than 1 on call to MA57C/CD. INFO(2) holds value of NRHS.\n",

                                 "LWORK too small on entry to MA57C/CD. INFO(2) holds minimum length\n"
                                 "required. A positive value of INFO(1) is associated with a warning\n"
                                 "message that will be output on unit ICNTL(2).\n"
                               };

  std::string ma57_wrn_msg[] = {
                                 "Operation successful.\n",

                                 "Index (in IRN or JCN) out of range on call to MA57A/AD or\n"
                                 "MA57D/DD. Action taken by subroutine is to ignore any such entries and\n"
                                 "continue. INFO(3) is set to the number of faulty entries. Details of the\n"
                                 "first ten are printed on unit ICNTL(2).\n",

                                 "Duplicate indices on call to MA57A/AD or MA57D/DD. Action taken by\n"
                                 "subroutine is to keep the duplicates and then to sum corresponding reals\n"
                                 "when MA57B/BD is called. INFO(4) is set to the number of faulty\n"
                                 "entries. Details of the first ten are printed on unit ICNTL(2).\n",

                                 "Both out-of-range indices and duplicates exist.\n",

                                 "Matrix is rank deficient on exit from MA57B/BD. In this case, a\n"
                                 "decomposition will still have been produced that will enable the\n"
                                 "subsequent solution of consistent equations. INFO(25) will be set to the\n"
                                 "rank of the factorized matrix.\n",

                                 "Pivots have different signs when factorizing a supposedly definite\n"
                                 "matrix (ICNTL(7) = 3) on call to MA57B/BD. INFO(26) is set to the number\n"
                                 "of sign changes.\n",

                                 "-",
                                 "-",

                                 "During error analysis the infinity norm of the computed solution was\n"
                                 "found to be zero.\n",

                                 "Insufficient real space to complete factorization when MA57B/BD called\n"
                                 "with ICNTL(8) != 0. User can copy real values to a longer array using\n"
                                 "MA57E/ED and recall MA57B/BD using this longer array to continue the\n"
                                 "factorization.\n",

                                 "Insufficient integer space to complete factorization when MA57B/BD\n"
                                 "called with ICNTL(8) != 0. User can copy integer values to a longer\n"
                                 "array using MA57E/ED and recall MA57B/BD using this longer array to\n"
                                 "continue the factorization.\n"
                               };

  Ma57TSolverInterface::Ma57TSolverInterface()
      :
      dim_(0),
      nonzeros_(0),
      initialized_(false),
      pivtol_changed_(false),
      refactorize_(false),

      wd_keep_(NULL),
      wd_iwork_(NULL),
      wd_fact_(NULL),
      wd_ifact_(NULL),
      a_(NULL)
  {
    DBG_START_METH("Ma57TSolverInterface::Ma57TSolverInterface()",
                   dbg_verbosity);
  }

  Ma57TSolverInterface::~Ma57TSolverInterface()
  {
    DBG_START_METH("Ma57TSolverInterface::~Ma57TSolverInterface()",
                   dbg_verbosity);
    delete [] a_;

    delete [] wd_fact_;
    delete [] wd_ifact_;

    delete [] wd_iwork_;
    delete [] wd_keep_;
  }

  void Ma57TSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption(
      "ma57_pivtol",
      "Pivot tolerance for the linear solver MA57.",
      0.0, true, 1.0, true, 1e-8,
      "A smaller number pivots for sparsity, a larger number pivots for "
      "stability. This option is only available if Ipopt has been compiled "
      "with MA57.");
    roptions->AddBoundedNumberOption(
      "ma57_pivtolmax",
      "Maximum pivot tolerance for the linear solver MA57.",
      0.0, true, 1.0, true, 1e-4,
      "Ipopt may increase pivtol as high as ma57_pivtolmax to get a more "
      "accurate solution to the linear system.  This option is only available "
      "if Ipopt has been compiled with MA57.");
    roptions->AddLowerBoundedNumberOption(
      "ma57_pre_alloc",
      "Safety factor for work space memory allocation for the linear solver MA57.",
      1., false, 3.,
      "If 1 is chosen, the suggested amount of work space is used.  However, "
      "choosing a larger number might avoid reallocation if the suggest values "
      "do not suffice.  This option is only available if Ipopt has been "
      "compiled with MA57.");
  }

  bool Ma57TSolverInterface::InitializeImpl(const OptionsList&  options,
      const std::string&    prefix)
  {
    // Obtain the options settings
    options.GetNumericValue("ma57_pivtol", pivtol_, prefix);
    if(options.GetNumericValue("ma57_pivtolmax", pivtolmax_, prefix)) {
      ASSERT_EXCEPTION(pivtolmax_>=pivtol_, OPTION_INVALID,
                       "Option \"pivtolmax\": This value must be between "
                       "pivtol and 1.");
    }
    else {
      pivtolmax_ = Max(pivtolmax_, pivtol_);
    }

    options.GetNumericValue("ma57_pre_alloc", ma57_pre_alloc_, prefix);

    // The following option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);
    DBG_ASSERT(!warm_start_same_structure_ && "warm_start_same_structure not yet implemented");


    /* Initialize. */
    F77_FUNC (ma57id, MA57ID) (wd_cntl_, wd_icntl_);

    /* Custom settings for MA57. */
    wd_icntl_[1-1] = 0;      /* Error stream */
    wd_icntl_[2-1] = 0;      /* Warning stream. */

    wd_icntl_[4-1] = 1;      /* Print statistics.  NOT Used. */
    wd_icntl_[5-1] = 0;      /* Print error. */

    // wd_icntl[6-1] = 0;       /* Pivoting order. */

    wd_cntl_[1-1]  = pivtol_;    /* Pivot threshold. */
    wd_icntl_[7-1] = 1;      /* Pivoting strategy. */

    // wd_icntl[8-1] = 0;       /* Retry factorization. */

    if (!warm_start_same_structure_) {
      dim_=0;
      nonzeros_=0;
      delete [] a_;
      delete [] wd_fact_;
      delete [] wd_ifact_;
      delete [] wd_iwork_;
      delete [] wd_keep_;
    }
    else {
      ASSERT_EXCEPTION(dim_>0 && nonzeros_>0, INVALID_WARMSTART,
                       "Ma57TSolverInterface called with warm_start_same_structure, "
                       "but the problem is solved for the first time.");
    }

    return true;
  }

  ESymSolverStatus
  Ma57TSolverInterface::MultiSolve(bool         new_matrix,
                                   const Index*     airn,
                                   const Index*     ajcn,
                                   Index        nrhs,
                                   double*      rhs_vals,
                                   bool         check_NegEVals,
                                   Index        numberOfNegEVals)
  {
    DBG_START_METH("Ma57TSolverInterface::MultiSolve",dbg_verbosity);

    // DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    // DBG_ASSERT(initialized_);
    // DBG_ASSERT(la_!=0);

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

  double* Ma57TSolverInterface::GetValuesArrayPtr()
  {
    DBG_START_METH("Ma57TSolverInterface::GetValuesArrayPtr",dbg_verbosity);
    DBG_ASSERT(initialized_);

    return a_;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus
  Ma57TSolverInterface::InitializeStructure(
    Index       dim,
    Index     nonzeros,
    const Index*  airn,
    const Index*  ajcn)
  {
    DBG_START_METH("Ma57TSolverInterface::InitializeStructure",dbg_verbosity);

    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    if (!warm_start_same_structure_) {
      dim_ = dim;
      nonzeros_ = nonzeros;
      // for MA57, a_ only has to be as long as the number of nonzero
      // elements
      delete [] a_;
      a_ = new double [nonzeros_];

      // Do the symbolic facotrization
      retval = SymbolicFactorization(airn, ajcn);
      if (retval != SYMSOLVER_SUCCESS ) {
        return retval;
      }
    }
    else {
      ASSERT_EXCEPTION(dim_==dim && nonzeros_==nonzeros, INVALID_WARMSTART,
                       "Ma57TSolverInterface called with warm_start_same_structure, "
                       "but the problem size has changed.");
    }

    initialized_ = true;

    return retval;
  }

  ESymSolverStatus
  Ma57TSolverInterface::SymbolicFactorization(const Index*  airn,
      const Index*  ajcn)
  {
    DBG_START_METH("Ma57TSolverInterface::SymbolicFactorization",dbg_verbosity);
    IpData().TimingStats().LinearSystemSymbolicFactorization().Start();

    ipfint n  = dim_;
    ipfint ne = nonzeros_;

    wd_lkeep_ = 5*n + ne + Max (n,ne) + 42;

    wd_cntl_[1-1]  = pivtol_;    /* Pivot threshold. */

    wd_iwork_ = new ipfint[5*n];
    wd_keep_  = new ipfint[wd_lkeep_];

    F77_FUNC (ma57ad, MA57AD)
    (&n, &ne, airn, ajcn, &wd_lkeep_, wd_keep_, wd_iwork_,
     wd_icntl_, wd_info_, wd_rinfo_);

    if (wd_info_[0] < 0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "*** Error from MA57AD *** INFO(0) = %d\n", wd_info_[0]);
    }

    wd_lfact_  = (ipfint)((Number)wd_info_[8] * ma57_pre_alloc_);
    wd_lifact_ = (ipfint)((Number)wd_info_[9] * ma57_pre_alloc_);

    // XXX MH:  Why is this necessary?  Is `::Factorization' called more
    // than once per object lifetime?  Where should allocation take
    // place, then?

    // AW: I moved this here now - my understanding is that wd_info[8]
    // and wd_info[9] are computed here, so we can just allocate the
    // amount of memory here.  I don't think there is any need to
    // reallocate it later for every factorization
    delete [] wd_fact_;
    delete [] wd_ifact_;

    wd_fact_  = new double[wd_lfact_];
    wd_ifact_ = new int[wd_lifact_];

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Suggested lfact  (*%e):  %d\n", ma57_pre_alloc_, wd_lfact_);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Suggested lifact (*%e):  %d\n", ma57_pre_alloc_, wd_lifact_);

    IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus
  Ma57TSolverInterface::Factorization(const Index*  airn,
                                      const Index*  ajcn,
                                      bool          check_NegEVals,
                                      Index         numberOfNegEVals)
  {
    DBG_START_METH("Ma57TSolverInterface::Factorization",dbg_verbosity);
    IpData().TimingStats().LinearSystemFactorization().Start();

    int fact_error = 1;

    wd_cntl_[1-1]  = pivtol_;    /* Pivot threshold. */

    ipfint n  = dim_;
    ipfint ne = nonzeros_;

    while (fact_error > 0) {
      F77_FUNC (ma57bd, MA57BD)
      (&n, &ne, a_, wd_fact_, &wd_lfact_, wd_ifact_, &wd_lifact_,
       &wd_lkeep_, wd_keep_, wd_iwork_,
       wd_icntl_, wd_cntl_, wd_info_, wd_rinfo_);

      negevals_ = wd_info_[24-1];   // Number of negative eigenvalues

      if (wd_info_[0] == 0) {
        fact_error = 0;
      }
      else if (wd_info_[0] == -3) {
        /* Failure due to insufficient REAL space on a call to MA57B/BD.
         * INFO(17) is set to a value that may suffice.  INFO(2) is set
         * to value of LFACT.  The user can allocate a larger array and
         * copy the contents of FACT into it using MA57E/ED, and recall
         * MA57B/BD.
         */
        double  *temp;
        ipfint  ic = 0;

        wd_lfact_ = wd_info_[16];
        temp = new double[wd_lfact_];

        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Reallocating lfact (%d)\n", wd_lfact_);

        ipfint idmy;
        F77_FUNC (ma57ed, MA57ED)
        (&n, &ic, wd_keep_,
         wd_fact_,  &wd_info_[1], temp, &wd_lfact_,
         wd_ifact_, &wd_info_[1], &idmy, &wd_lfact_,
         wd_info_);

        delete [] wd_fact_;
        wd_fact_ = temp;
      }
      else if (wd_info_[0] == -4) {
        /* Failure due to insufficient INTEGER space on a call to
         * MA57B/BD.  INFO(18) is set to a value that may suffice.
         * INFO(2) is set to value of LIFACT.  The user can allocate a
         * larger array and copy the contents of IFACT into it using
         * MA57E/ED, and recall MA57B/BD.
         */

        ipfint  *temp;
        ipfint   ic = 1;

        wd_lifact_ = wd_info_[17];
        temp = new ipfint[wd_lifact_];

        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "Reallocating lifact (%d)\n", wd_lifact_);

        double ddmy;
        F77_FUNC (ma57ed, MA57ED)
        (&n, &ic, wd_keep_,
         wd_fact_,  &wd_info_[1], &ddmy, &wd_lifact_,
         wd_ifact_, &wd_info_[1], temp, &wd_lifact_,
         wd_info_);

        delete [] wd_ifact_;
        wd_ifact_ = temp;
      }
      else if (wd_info_[0] < 0) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error in MA57BD:  %d\n", wd_info_[0]);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "MA57 Error message: %s\n",
                       ma57_err_msg[-wd_info_[1-1]].c_str());
        return SYMSOLVER_FATAL_ERROR;
      }
      // Check if the system is singular.
      else if (wd_info_[0] == 4) {
        IpData().TimingStats().LinearSystemFactorization().End();
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "System singular, rank = %d\n", wd_info_[25-1]);
        return SYMSOLVER_SINGULAR;
      }
      else if (wd_info_[0] > 0) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Warning in MA57BD:  %d\n", wd_info_[0]);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "MA57 Warning message: %s\n",
                       ma57_wrn_msg[wd_info_[1-1]].c_str());
        // For now, abort the process so that we don't miss any problems
        return SYMSOLVER_FATAL_ERROR;
      }
    }

    // double peak_mem = 1.0e-3 * (wd_lfact*8.0 + wd_lifact*4.0 + wd_lkeep*4.0);


    // Check whether the number of negative eigenvalues matches the
    // requested count.
    IpData().TimingStats().LinearSystemFactorization().End();
    if (check_NegEVals && (numberOfNegEVals!=negevals_)) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma57TSolverInterface::Factorization: "
                     "negevals_ = %d, but numberOfNegEVals = %d\n",
                     negevals_, numberOfNegEVals);
      return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
  }

  ESymSolverStatus Ma57TSolverInterface::Backsolve(
    Index     nrhs,
    double    *rhs_vals)
  {
    DBG_START_METH("Ma27TSolverInterface::Backsolve",dbg_verbosity);
    IpData().TimingStats().LinearSystemBackSolve().Start();

    ipfint  n      = dim_;
    ipfint  job    = 1;

    ipfint  nrhs_X = nrhs;
    ipfint  lrhs   = n;

    ipfint  lwork;
    double* work;

    lwork = n * nrhs;
    work = new double[lwork];

    // For each right hand side, call MA57CD
    // XXX MH: MA57 can do several RHSs; just do one solve...
    // AW: Ok is the following correct?
    if (DBG_VERBOSITY()>=2) {
      for(Index irhs=0; irhs<nrhs; irhs++) {
        for (Index i=0; i<dim_; i++) {
          DBG_PRINT((2, "rhs[%2d,%5d] = %23.15e\n", irhs, i, rhs_vals[irhs*dim_+i]));
        }
      }
    }

    F77_FUNC (ma57cd, MA57CD)
    (&job, &n, wd_fact_, &wd_lfact_, wd_ifact_, &wd_lifact_,
     &nrhs_X, rhs_vals, &lrhs,
     work, &lwork, wd_iwork_,
     wd_icntl_, wd_info_);

    if (wd_info_[0] != 0)
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in MA57CD:  %d.\n", wd_info_[0]);

    if (DBG_VERBOSITY()>=2) {
      for(Index irhs=0; irhs<nrhs; irhs++) {
        for (Index i=0; i<dim_; i++) {
          DBG_PRINT((2, "sol[%2d,%5d] = %23.15e\n", irhs, i, rhs_vals[irhs*dim_+i]));
        }
      }
    }

    delete [] work;

    IpData().TimingStats().LinearSystemBackSolve().End();
    return SYMSOLVER_SUCCESS;
  }

  Index Ma57TSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("Ma57TSolverInterface::NumberOfNegEVals",dbg_verbosity);
    DBG_ASSERT(ProvidesInertia());
    DBG_ASSERT(initialized_);
    return negevals_;
  }

  bool Ma57TSolverInterface::IncreaseQuality()
  {
    DBG_START_METH("Ma57TSolverInterface::IncreaseQuality",dbg_verbosity);
    if (pivtol_ == pivtolmax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Indreasing pivot tolerance for MA57 from %7.2e ",
                   pivtol_);
    pivtol_ = Min(pivtolmax_, pow(pivtol_,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   pivtol_);
    return true;
  }

} // namespace Ipopt
