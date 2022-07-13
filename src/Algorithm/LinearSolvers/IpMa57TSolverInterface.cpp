// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Michael Hagemann               Univ of Basel 2005-10-28
//               original version (based on MA27TSolverInterface.cpp)

#include "IpoptConfig.h"
#include "IpMa57TSolverInterface.hpp"

#include <cmath>
#include <iostream>

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

#if (defined(COINHSL_HAS_MA57) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA57S) && defined(IPOPT_SINGLE))
#ifdef IPOPT_SINGLE
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name,NAME)
#else
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name ## d,NAME ## D)
#endif

/** MA57 functions from HSL library (symbols resolved at linktime) */
extern "C"
{
   IPOPT_DECL_MA57A(IPOPT_HSL_FUNCP(ma57a, MA57A));
   IPOPT_DECL_MA57B(IPOPT_HSL_FUNCP(ma57b, MA57B));
   IPOPT_DECL_MA57C(IPOPT_HSL_FUNCP(ma57c, MA57C));
   IPOPT_DECL_MA57E(IPOPT_HSL_FUNCP(ma57e, MA57E));
   IPOPT_DECL_MA57I(IPOPT_HSL_FUNCP(ma57i, MA57I));
}
#else
#ifdef IPOPT_SINGLE
#define HSLFUNCNAMESUFFIX ""
#else
#define HSLFUNCNAMESUFFIX "d"
#endif
#endif

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

/** pointer to MA57 function that can be set via Ma57TSolverInterface::SetFunctions() */
static IPOPT_DECL_MA57A(*user_ma57a) = NULL;
static IPOPT_DECL_MA57B(*user_ma57b) = NULL;
static IPOPT_DECL_MA57C(*user_ma57c) = NULL;
static IPOPT_DECL_MA57E(*user_ma57e) = NULL;
static IPOPT_DECL_MA57I(*user_ma57i) = NULL;

const char* ma57_err_msg[] =
{
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

const char* ma57_wrn_msg[] =
{
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

   "-", "-",

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

Ma57TSolverInterface::Ma57TSolverInterface(
   SmartPtr<LibraryLoader> hslloader_
)  : hslloader(hslloader_),
   ma57a(NULL),
   ma57b(NULL),
   ma57c(NULL),
   ma57e(NULL),
   ma57i(NULL),
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
   DBG_START_METH("Ma57TSolverInterface::Ma57TSolverInterface()", dbg_verbosity);
}

Ma57TSolverInterface::~Ma57TSolverInterface()
{
   DBG_START_METH("Ma57TSolverInterface::~Ma57TSolverInterface()", dbg_verbosity);
   delete[] a_;

   delete[] wd_fact_;
   delete[] wd_ifact_;

   delete[] wd_iwork_;
   delete[] wd_keep_;
}

void Ma57TSolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddLowerBoundedIntegerOption(
      "ma57_print_level",
      "Debug printing level for the linear solver MA57",
      0, 0,
      "0: no printing; 1: Error messages only; 2: Error and warning messages; 3: Error and warning messages and terse monitoring; >=4: All information.");
   roptions->AddBoundedNumberOption(
      "ma57_pivtol",
      "Pivot tolerance for the linear solver MA57.",
      0.0, true,
      1.0, true,
      1e-8,
      "A smaller number pivots for sparsity, a larger number pivots for stability.");
   roptions->AddBoundedNumberOption(
      "ma57_pivtolmax",
      "Maximum pivot tolerance for the linear solver MA57.",
      0.0, true,
      1.0, true,
      1e-4,
      "Ipopt may increase pivtol as high as ma57_pivtolmax to get a more accurate solution to the linear system.");
   roptions->AddLowerBoundedNumberOption(
      "ma57_pre_alloc",
      "Safety factor for work space memory allocation for the linear solver MA57.",
      1., false,
      1.05,
      "If 1 is chosen, the suggested amount of work space is used. "
      "However, choosing a larger number might avoid reallocation if the suggest values do not suffice.");
   roptions->AddBoundedIntegerOption(
      "ma57_pivot_order",
      "Controls pivot order in MA57",
      0, 5,
#ifdef FUNNY_MA57_FINT
      2, // Matlab's MA57 can crash if you try to use Metis
#else
      5,
#endif
      "This is ICNTL(6) in MA57.");
   roptions->AddBoolOption(
      "ma57_automatic_scaling",
      "Controls whether to enable automatic scaling in MA57",
      false,
      "For higher reliability of the MA57 solver, you may want to set this option to yes. "
      "This is ICNTL(15) in MA57.");

   // CET: 04-29-2010
   roptions->AddLowerBoundedIntegerOption(
      "ma57_block_size",
      "Controls block size used by Level 3 BLAS in MA57BD",
      1,
      16,
      "This is ICNTL(11) in MA57.");

   roptions->AddLowerBoundedIntegerOption(
      "ma57_node_amalgamation",
      "Node amalgamation parameter",
      1,
      16,
      "This is ICNTL(12) in MA57.");

   roptions->AddBoundedIntegerOption(
      "ma57_small_pivot_flag",
      "Handling of small pivots",
      0, 1,
      0,
      "If set to 1, then when small entries defined by CNTL(2) are detected they are removed and "
      "the corresponding pivots placed at the end of the factorization. "
      "This can be particularly efficient if the matrix is highly rank deficient. "
      "This is ICNTL(16) in MA57.");
   // CET 04-29-2010
}

/// set MA57 functions to use for every instantiation of this class
void Ma57TSolverInterface::SetFunctions(
   IPOPT_DECL_MA57A(*ma57a),
   IPOPT_DECL_MA57B(*ma57b),
   IPOPT_DECL_MA57C(*ma57c),
   IPOPT_DECL_MA57E(*ma57e),
   IPOPT_DECL_MA57I(*ma57i)
)
{
   DBG_ASSERT(ma57a != NULL);
   DBG_ASSERT(ma57b != NULL);
   DBG_ASSERT(ma57c != NULL);
   DBG_ASSERT(ma57e != NULL);
   DBG_ASSERT(ma57i != NULL);

   user_ma57a = ma57a;
   user_ma57b = ma57b;
   user_ma57c = ma57c;
   user_ma57e = ma57e;
   user_ma57i = ma57i;
}

bool Ma57TSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( user_ma57a != NULL )
   {
      // someone set MA57 functions via setFunctions - prefer these
      ma57a = user_ma57a;
      ma57b = user_ma57b;
      ma57c = user_ma57c;
      ma57e = user_ma57e;
      ma57i = user_ma57i;
   }
   else
   {
#if (defined(COINHSL_HAS_MA57) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA57S) && defined(IPOPT_SINGLE))
      // use HSL functions that should be available in linked HSL library
      ma57a = &::IPOPT_HSL_FUNCP(ma57a, MA57A);
      ma57b = &::IPOPT_HSL_FUNCP(ma57b, MA57B);
      ma57c = &::IPOPT_HSL_FUNCP(ma57c, MA57C);
      ma57e = &::IPOPT_HSL_FUNCP(ma57e, MA57E);
      ma57i = &::IPOPT_HSL_FUNCP(ma57i, MA57I);
#else
      // try to load HSL functions from a shared library at runtime
      DBG_ASSERT(IsValid(hslloader));

      ma57a = (IPOPT_DECL_MA57A(*))hslloader->loadSymbol("ma57a" HSLFUNCNAMESUFFIX);
      ma57b = (IPOPT_DECL_MA57B(*))hslloader->loadSymbol("ma57b" HSLFUNCNAMESUFFIX);
      ma57c = (IPOPT_DECL_MA57C(*))hslloader->loadSymbol("ma57c" HSLFUNCNAMESUFFIX);
      ma57e = (IPOPT_DECL_MA57E(*))hslloader->loadSymbol("ma57e" HSLFUNCNAMESUFFIX);
      ma57i = (IPOPT_DECL_MA57I(*))hslloader->loadSymbol("ma57i" HSLFUNCNAMESUFFIX);
#endif
   }

   DBG_ASSERT(ma57a != NULL);
   DBG_ASSERT(ma57b != NULL);
   DBG_ASSERT(ma57c != NULL);
   DBG_ASSERT(ma57e != NULL);
   DBG_ASSERT(ma57i != NULL);

   // Obtain the options settings
   Index print_level;
   options.GetIntegerValue("ma57_print_level", print_level, prefix);

   options.GetNumericValue("ma57_pivtol", pivtol_, prefix);
   if( options.GetNumericValue("ma57_pivtolmax", pivtolmax_, prefix) )
   {
      ASSERT_EXCEPTION(pivtolmax_ >= pivtol_, OPTION_INVALID, "Option \"pivtolmax\": This value must be between "
                       "pivtol and 1.");
   }
   else if( pivtol_ > pivtolmax_ )
   {
      pivtolmax_ = pivtol_;
   }

   options.GetNumericValue("ma57_pre_alloc", ma57_pre_alloc_, prefix);
   Index ma57_pivot_order;
   options.GetIntegerValue("ma57_pivot_order", ma57_pivot_order, prefix);

   // The following option is registered by OrigIpoptNLP
   options.GetBoolValue("warm_start_same_structure", warm_start_same_structure_, prefix);
   DBG_ASSERT(!warm_start_same_structure_ && "warm_start_same_structure not yet implemented");

   bool ma57_automatic_scaling;
   options.GetBoolValue("ma57_automatic_scaling", ma57_automatic_scaling, prefix);

   // CET 04-29-2010
   Index ma57_block_size;
   options.GetIntegerValue("ma57_block_size", ma57_block_size, prefix);

   Index ma57_node_amalgamation;
   options.GetIntegerValue("ma57_node_amalgamation", ma57_node_amalgamation, prefix);

   Index ma57_small_pivot_flag;
   options.GetIntegerValue("ma57_small_pivot_flag", ma57_small_pivot_flag, prefix);
   // CET 04-29-2010

   /* Initialize. */
   ma57i(wd_cntl_, wd_icntl_);
   /* Custom settings for MA57. */
   wd_icntl_[0] = 0; /* Error stream */
   wd_icntl_[1] = 0; /* Warning stream. */

   wd_icntl_[3] = 1; /* Print statistics.  NOT Used. */
   wd_icntl_[4] = print_level; /* Print level. */

   wd_icntl_[5] = ma57_pivot_order; /* Pivoting order. */

   wd_cntl_[0] = pivtol_; /* Pivot threshold. */
   wd_icntl_[6] = 1; /* Pivoting strategy. */

   // CET: Added 04-29-2010 at suggestion of Jonathan Hogg of HSL
   wd_icntl_[10] = ma57_block_size; /* Block size used by Level 3 BLAS in MA57BD - should be a multiple of 8.  Default is 16. */
   wd_icntl_[11] = ma57_node_amalgamation; /* Two nodes of the assembly tree are merged only if both involve less than ICNTL(12) eliminations.  Default is 16. */
   // CET: 04-29-2010

   if( ma57_automatic_scaling )
   {
      wd_icntl_[14] = 1;
   }
   else
   {
      wd_icntl_[14] = 0;
   }

   // CET: Added 04-29-2010 at suggestion of Jonathan Hogg of HSL
   wd_icntl_[15] = ma57_small_pivot_flag; /* If set to 1, small entries are removed and corresponding pivots are placed at the end of factorization.  May be useful for highly rank deficient matrices.  Default is 0. */
   // CET: 04-29-2010

   // wd_icntl[8-1] = 0;       /* Retry factorization. */

   if( !warm_start_same_structure_ )
   {
      dim_ = 0;
      nonzeros_ = 0;
      delete[] a_;
      a_ = NULL;
      delete[] wd_fact_;
      wd_fact_ = NULL;
      delete[] wd_ifact_;
      wd_ifact_ = NULL;
      delete[] wd_iwork_;
      wd_iwork_ = NULL;
      delete[] wd_keep_;
      wd_keep_ = NULL;
   }
   else
   {
      ASSERT_EXCEPTION(dim_ > 0 && nonzeros_ > 0, INVALID_WARMSTART,
                       "Ma57TSolverInterface called with warm_start_same_structure, "
                       "but the problem is solved for the first time.");
   }

   return true;
}

ESymSolverStatus Ma57TSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* airn,
   const Index* ajcn,
   Index        nrhs,
   Number*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("Ma57TSolverInterface::MultiSolve", dbg_verbosity);

   // DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   // DBG_ASSERT(initialized_);
   // DBG_ASSERT(la_!=0);

   if( pivtol_changed_ )
   {
      DBG_PRINT((1, "Pivot tolerance has changed.\n"));
      pivtol_changed_ = false;
      // If the pivot tolerance has been changed but the matrix is not
      // new, we have to request the values for the matrix again to do
      // the factorization again.
      if( !new_matrix )
      {
         DBG_PRINT((1, "Ask caller to call again.\n"));
         refactorize_ = true;
         return SYMSOLVER_CALL_AGAIN;
      }
   }

   // check if a factorization has to be done
   DBG_PRINT((1, "new_matrix = %d\n", new_matrix));
   if( new_matrix || refactorize_ )
   {
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(airn, ajcn, check_NegEVals, numberOfNegEVals);
      if( retval != SYMSOLVER_SUCCESS )
      {
         DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
         return retval;  // Matrix singular or error occurred
      }
      refactorize_ = false;
   }

   // do the backsolve
   return Backsolve(nrhs, rhs_vals);
}

Number* Ma57TSolverInterface::GetValuesArrayPtr()
{
   DBG_START_METH("Ma57TSolverInterface::GetValuesArrayPtr", dbg_verbosity);
   DBG_ASSERT(initialized_);

   return a_;
}

/** Initialize the local copy of the positions of the nonzero elements */
ESymSolverStatus Ma57TSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* airn,
   const Index* ajcn
)
{
   DBG_START_METH("Ma57TSolverInterface::InitializeStructure", dbg_verbosity);

   ESymSolverStatus retval = SYMSOLVER_SUCCESS;
   if( !warm_start_same_structure_ )
   {
      dim_ = dim;
      nonzeros_ = nonzeros;
      // for MA57, a_ only has to be as long as the number of nonzero
      // elements
      delete[] a_;
      a_ = NULL;
      a_ = new Number[nonzeros_];

      // Do the symbolic factorization
      retval = SymbolicFactorization(airn, ajcn);
      if( retval != SYMSOLVER_SUCCESS )
      {
         return retval;
      }
   }
   else
   {
      ASSERT_EXCEPTION(dim_ == dim && nonzeros_ == nonzeros, INVALID_WARMSTART,
                       "Ma57TSolverInterface called with warm_start_same_structure, "
                       "but the problem size has changed.");
   }

   initialized_ = true;

   return retval;
}

ESymSolverStatus Ma57TSolverInterface::SymbolicFactorization(
   const Index* airn,
   const Index* ajcn
)
{
   DBG_START_METH("Ma57TSolverInterface::SymbolicFactorization", dbg_verbosity);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
   }

   ma57int n = dim_;
   ma57int ne = nonzeros_;

   wd_lkeep_ = 5 * n + ne + (n > ne ? n : ne) + 42;

   wd_cntl_[0] = pivtol_; /* Pivot threshold. */

   wd_iwork_ = new ma57int[5 * n];
   wd_keep_ = new ma57int[wd_lkeep_];
   // Initialize to 0 as otherwise MA57EX can sometimes fail
   for( ma57int k = 0; k < wd_lkeep_; k++ )
   {
      wd_keep_[k] = 0;
   }

   // copy-cast airn and ajcn into ma57int arrays
   ma57int* airn_ma57int;
   ma57int* ajcn_ma57int;
   if( sizeof(ma57int) != sizeof(Index) )
   {
      airn_ma57int = new ma57int[ne];
      ajcn_ma57int = new ma57int[ne];
      for( ma57int k = 0; k < ne; ++k )
      {
         airn_ma57int[k] = (ma57int) airn[k];
         ajcn_ma57int[k] = (ma57int) ajcn[k];
      }
   }
   else
   {
      airn_ma57int = (ma57int*) (void*) const_cast<Index*>(airn);
      ajcn_ma57int = (ma57int*) (void*) const_cast<Index*>(ajcn);
   }
   ma57a(&n, &ne, airn_ma57int, ajcn_ma57int, &wd_lkeep_, wd_keep_, wd_iwork_, wd_icntl_, wd_info_, wd_rinfo_);
   // free copy-casted ma57int arrays, no longer needed
   if( sizeof(ma57int) != sizeof(Index) )
   {
      delete[] airn_ma57int;
      delete[] ajcn_ma57int;
   }

   if( wd_info_[0] < 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "*** Error from MA57AD *** INFO(0) = %" IPOPT_INDEX_FORMAT "\n", wd_info_[0]);
   }

   wd_lfact_ = 0;
   wd_lifact_ = 0;
   ComputeMemIncrease(wd_lfact_, (Number)wd_info_[8] * ma57_pre_alloc_, 0, "double working space for MA57");
   ComputeMemIncrease(wd_lifact_, (Number)wd_info_[9] * ma57_pre_alloc_, 0, "integer working space for MA57");

   // XXX MH:  Why is this necessary?  Is `::Factorization' called more
   // than once per object lifetime?  Where should allocation take
   // place, then?

   // AW: I moved this here now - my understanding is that wd_info[8]
   // and wd_info[9] are computed here, so we can just allocate the
   // amount of memory here.  I don't think there is any need to
   // reallocate it later for every factorization
   delete[] wd_fact_;
   wd_fact_ = NULL;
   delete[] wd_ifact_;
   wd_ifact_ = NULL;

   wd_fact_ = new Number[wd_lfact_];
   wd_ifact_ = new ma57int[wd_lifact_];

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Suggested lfact  (*%e):  %" IPOPT_INDEX_FORMAT "\n", ma57_pre_alloc_, wd_lfact_);
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Suggested lifact (*%e):  %" IPOPT_INDEX_FORMAT "\n", ma57_pre_alloc_, wd_lifact_);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
   }
   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus Ma57TSolverInterface::Factorization(
   const Index* /*airn*/,
   const Index* /*ajcn*/,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("Ma57TSolverInterface::Factorization", dbg_verbosity);
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().Start();
   }

   bool fact_error = true;

   wd_cntl_[0] = pivtol_; /* Pivot threshold. */

   ma57int n = dim_;
   ma57int ne = nonzeros_;

   while( fact_error )
   {
      ma57b(&n, &ne, a_, wd_fact_, &wd_lfact_, wd_ifact_, &wd_lifact_, &wd_lkeep_, wd_keep_,
            wd_iwork_, wd_icntl_, wd_cntl_, wd_info_, wd_rinfo_);
      negevals_ = (Index) wd_info_[24 - 1]; // Number of negative eigenvalues

      if( wd_info_[0] == 0 )
      {
         fact_error = false;
      }
      else if( wd_info_[0] == -3 )
      {
         /* Failure due to insufficient REAL space on a call to MA57B/BD.
          * INFO(17) is set to a value that may suffice.  INFO(2) is set
          * to value of LFACT.  The user can allocate a larger array and
          * copy the contents of FACT into it using MA57E/ED, and recall
          * MA57B/BD.
          */
         Number* temp;
         ma57int ic = 0;

         ComputeMemIncrease(wd_lfact_, (Number)wd_info_[16] * ma57_pre_alloc_, 0, "double working space for MA57");
         Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                        "Reallocating memory for MA57: lfact (%" IPOPT_INDEX_FORMAT ")\n", wd_lfact_);

         // I removed this, because the call to new below should already check this, too, and would throw a std::bad_alloc (which we catch) (or std::bad_array_new_length if C++11)
         // if( (size_t) wd_lfact_ > std::numeric_limits<size_t>::max() / sizeof(Number) )
         // {
         //   Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
         //                  "Cannot allocate memory of size %" IPOPT_INDEX_FORMAT " exceeding SIZE_MAX = %zd\n", wd_lfact_, std::numeric_limits<size_t>::max());
         //   return SYMSOLVER_FATAL_ERROR;
         // }

         temp = new Number[wd_lfact_];

         ma57int idmy;
         ma57e(&n, &ic, wd_keep_, wd_fact_, &wd_info_[1], temp, &wd_lfact_, wd_ifact_, &wd_info_[1],
               &idmy, &wd_lfact_, wd_info_);
         delete[] wd_fact_;
         wd_fact_ = temp;
      }
      else if( wd_info_[0] == -4 )
      {
         /* Failure due to insufficient INTEGER space on a call to
          * MA57B/BD.  INFO(18) is set to a value that may suffice.
          * INFO(2) is set to value of LIFACT.  The user can allocate a
          * larger array and copy the contents of IFACT into it using
          * MA57E/ED, and recall MA57B/BD.
          */

         ma57int* temp;
         ma57int ic = 1;

         ComputeMemIncrease(wd_lifact_, (Number)wd_info_[17] * ma57_pre_alloc_, 0, "integer working space for MA57");
         temp = new ma57int[wd_lifact_];

         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Reallocating lifact (%" IPOPT_INDEX_FORMAT ")\n", wd_lifact_);

         Number ddmy;
         ma57e(&n, &ic, wd_keep_, wd_fact_, &wd_info_[1], &ddmy, &wd_lifact_, wd_ifact_,
               &wd_info_[1], temp, &wd_lifact_, wd_info_);
         delete[] wd_ifact_;
         wd_ifact_ = temp;
      }
      else if( wd_info_[0] < 0 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error in MA57BD:  %" IPOPT_INDEX_FORMAT "\n", wd_info_[0]);
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "MA57 Error message: %s\n", ma57_err_msg[-wd_info_[0]]);
         return SYMSOLVER_FATAL_ERROR;
      }
      // Check if the system is singular.
      else if( wd_info_[0] == 4 )
      {
         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemFactorization().End();
         }
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "System singular, rank = %" IPOPT_INDEX_FORMAT "\n", wd_info_[24]);
         return SYMSOLVER_SINGULAR;
      }
      else if( wd_info_[0] > 0 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Warning in MA57BD:  %" IPOPT_INDEX_FORMAT "\n", wd_info_[0]);
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "MA57 Warning message: %s\n", ma57_wrn_msg[wd_info_[0]]);
         // For now, abort the process so that we don't miss any problems
         return SYMSOLVER_FATAL_ERROR;
      }
   }

   Number peak_mem = 1.0e-3 * ((Number) wd_lfact_ * sizeof(Number) + (Number) wd_lifact_ * sizeof(ma57int) + (Number) wd_lkeep_ * sizeof(ma57int));
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "MA57 peak memory use: %zdKB\n", (size_t) (peak_mem));

   // Check whether the number of negative eigenvalues matches the
   // requested count.
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().End();
   }
   if( check_NegEVals && (numberOfNegEVals != negevals_) )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma57TSolverInterface::Factorization: negevals_ = %" IPOPT_INDEX_FORMAT ", but numberOfNegEVals = %" IPOPT_INDEX_FORMAT "\n", negevals_, numberOfNegEVals);
      return SYMSOLVER_WRONG_INERTIA;
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus Ma57TSolverInterface::Backsolve(
   Index   nrhs,
   Number* rhs_vals
)
{
   DBG_START_METH("Ma57TSolverInterface::Backsolve", dbg_verbosity);
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().Start();
   }

   ma57int n = dim_;
   ma57int job = 1;

   ma57int nrhs_X = nrhs;
   ma57int lrhs = n;

   ma57int lwork;
   Number* work;

   lwork = n * nrhs;
   work = new Number[lwork];

   // For each right hand side, call MA57CD
   // XXX MH: MA57 can do several RHSs; just do one solve...
   // AW: Ok is the following correct?
   if( DBG_VERBOSITY() >= 2 )
   {
      for( Index irhs = 0; irhs < nrhs; irhs++ )
      {
         for( Index i = 0; i < dim_; i++ )
         {
            DBG_PRINT((2, "rhs[%2d,%5d] = %23.15e\n", irhs, i, rhs_vals[irhs * dim_ + i]));
         }
      }
   }
   ma57c(&job, &n, wd_fact_, &wd_lfact_, wd_ifact_, &wd_lifact_, &nrhs_X, rhs_vals, &lrhs, work,
         &lwork, wd_iwork_, wd_icntl_, wd_info_);
   if( wd_info_[0] != 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in MA57CD:  %" IPOPT_INDEX_FORMAT ".\n", wd_info_[0]);
   }

   if( DBG_VERBOSITY() >= 2 )
   {
      for( Index irhs = 0; irhs < nrhs; irhs++ )
      {
         for( Index i = 0; i < dim_; i++ )
         {
            DBG_PRINT((2, "sol[%2d,%5d] = %23.15e\n", irhs, i, rhs_vals[irhs * dim_ + i]));
         }
      }
   }

   delete[] work;

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().End();
   }
   return SYMSOLVER_SUCCESS;
}

Index Ma57TSolverInterface::NumberOfNegEVals() const
{
   DBG_START_METH("Ma57TSolverInterface::NumberOfNegEVals", dbg_verbosity);
   DBG_ASSERT(ProvidesInertia());
   DBG_ASSERT(initialized_);
   return negevals_;
}

bool Ma57TSolverInterface::IncreaseQuality()
{
   DBG_START_METH("Ma57TSolverInterface::IncreaseQuality", dbg_verbosity);
   if( pivtol_ == pivtolmax_ )
   {
      return false;
   }
   pivtol_changed_ = true;

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Increasing pivot tolerance for MA57 from %7.2e ", pivtol_);
   pivtol_ = Min(pivtolmax_, std::pow(pivtol_, Number(0.75)));
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "to %7.2e.\n", pivtol_);
   return true;
}

} // namespace Ipopt
