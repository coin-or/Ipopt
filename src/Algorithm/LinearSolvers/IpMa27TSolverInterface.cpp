// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#include "IpoptConfig.h"
#include "IpMa27TSolverInterface.hpp"

#include <cmath>

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

#if (defined(COINHSL_HAS_MA27) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA27S) && defined(IPOPT_SINGLE))
#ifdef IPOPT_SINGLE
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name,NAME)
#else
#define IPOPT_HSL_FUNCP(name,NAME) IPOPT_HSL_FUNC(name ## d,NAME ## D)
#endif

/** MA27 functions from HSL library (symbols resolved at linktime) */
extern "C"
{
   IPOPT_DECL_MA27A(IPOPT_HSL_FUNCP(ma27a, MA27A));
   IPOPT_DECL_MA27B(IPOPT_HSL_FUNCP(ma27b, MA27B));
   IPOPT_DECL_MA27C(IPOPT_HSL_FUNCP(ma27c, MA27C));
   IPOPT_DECL_MA27I(IPOPT_HSL_FUNCP(ma27i, MA27I));
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

/** pointer to MA27 function that can be set via Ma27TSolverInterface::SetFunctions() */
static IPOPT_DECL_MA27A(*user_ma27a) = NULL;
static IPOPT_DECL_MA27B(*user_ma27b) = NULL;
static IPOPT_DECL_MA27C(*user_ma27c) = NULL;
static IPOPT_DECL_MA27I(*user_ma27i) = NULL;

Ma27TSolverInterface::Ma27TSolverInterface(
   SmartPtr<LibraryLoader> hslloader_
)  : hslloader(hslloader_),
   ma27a(NULL),
   ma27b(NULL),
   ma27c(NULL),
   ma27i(NULL),
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
   DBG_START_METH("Ma27TSolverInterface::Ma27TSolverInterface()", dbg_verbosity);
}

Ma27TSolverInterface::~Ma27TSolverInterface()
{
   DBG_START_METH("Ma27TSolverInterface::~Ma27TSolverInterface()",
                  dbg_verbosity);
   delete[] iw_;
   delete[] ikeep_;
   delete[] a_;
}

void Ma27TSolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddBoundedIntegerOption(
      "ma27_print_level",
      "Debug printing level for the linear solver MA27",
      0, 4, 0,
      "0: no printing; 1: Error messages only; 2: Error and warning messages; 3: Error and warning messages and terse monitoring; 4: All information.");
   roptions->AddBoundedNumberOption(
      "ma27_pivtol",
      "Pivot tolerance for the linear solver MA27.",
      0.0, true,
      1.0, true,
      1e-8,
      "A smaller number pivots for sparsity, a larger number pivots for stability.");
   roptions->AddBoundedNumberOption(
      "ma27_pivtolmax",
      "Maximum pivot tolerance for the linear solver MA27.",
      0.0, true,
      1.0, true,
      1e-4,
      "Ipopt may increase pivtol as high as ma27_pivtolmax to get a more accurate solution to the linear system.");
   roptions->AddLowerBoundedNumberOption(
      "ma27_liw_init_factor",
      "Integer workspace memory for MA27.",
      1.0, false,
      5.0,
      "The initial integer workspace memory = liw_init_factor * memory required by unfactored system. "
      "Ipopt will increase the workspace size by ma27_meminc_factor if required.");
   roptions->AddLowerBoundedNumberOption(
      "ma27_la_init_factor",
      "Real workspace memory for MA27.",
      1.0, false,
      5.0,
      "The initial real workspace memory = la_init_factor * memory required by unfactored system. "
      "Ipopt will increase the workspace size by ma27_meminc_factor if required.");
   roptions->AddLowerBoundedNumberOption(
      "ma27_meminc_factor",
      "Increment factor for workspace size for MA27.",
      1.0, false,
      2.0,
      "If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.");
   roptions->AddBoolOption(
      "ma27_skip_inertia_check",
      "Whether to always pretend that inertia is correct.",
      false,
      "Setting this option to \"yes\" essentially disables inertia check. "
      "This option makes the algorithm non-robust and easily fail, but it might give some insight into the necessity of inertia control.",
      true);
   roptions->AddBoolOption(
      "ma27_ignore_singularity",
      "Whether to use MA27's ability to solve a linear system even if the matrix is singular.",
      false,
      "Setting this option to \"yes\" means that Ipopt will call MA27 to compute solutions for right hand sides, "
      "even if MA27 has detected that the matrix is singular (but is still able to solve the linear system). "
      "In some cases this might be better than using Ipopt's heuristic of small perturbation of the lower diagonal of the KKT matrix.",
      true);
}

void Ma27TSolverInterface::SetFunctions(
   IPOPT_DECL_MA27A(*ma27a),
   IPOPT_DECL_MA27B(*ma27b),
   IPOPT_DECL_MA27C(*ma27c),
   IPOPT_DECL_MA27I(*ma27i)
)
{
   DBG_ASSERT(ma27a != NULL);
   DBG_ASSERT(ma27b != NULL);
   DBG_ASSERT(ma27c != NULL);
   DBG_ASSERT(ma27i != NULL);

   user_ma27a = ma27a;
   user_ma27b = ma27b;
   user_ma27c = ma27c;
   user_ma27i = ma27i;
}

bool Ma27TSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( user_ma27a != NULL )
   {
      // someone set MA27 functions via setFunctions - prefer these
      ma27a = user_ma27a;
      ma27b = user_ma27b;
      ma27c = user_ma27c;
      ma27i = user_ma27i;
   }
   else
   {
#if (defined(COINHSL_HAS_MA27) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA27S) && defined(IPOPT_SINGLE))
      // use HSL functions that should be available in linked HSL library
      ma27a = &::IPOPT_HSL_FUNCP(ma27a, MA27A);
      ma27b = &::IPOPT_HSL_FUNCP(ma27b, MA27B);
      ma27c = &::IPOPT_HSL_FUNCP(ma27c, MA27C);
      ma27i = &::IPOPT_HSL_FUNCP(ma27i, MA27I);
#else
      DBG_ASSERT(IsValid(hslloader));

      ma27a = (IPOPT_DECL_MA27A(*))hslloader->loadSymbol("ma27a" HSLFUNCNAMESUFFIX);
      ma27b = (IPOPT_DECL_MA27B(*))hslloader->loadSymbol("ma27b" HSLFUNCNAMESUFFIX);
      ma27c = (IPOPT_DECL_MA27C(*))hslloader->loadSymbol("ma27c" HSLFUNCNAMESUFFIX);
      ma27i = (IPOPT_DECL_MA27I(*))hslloader->loadSymbol("ma27i" HSLFUNCNAMESUFFIX);
#endif
   }

   DBG_ASSERT(ma27a != NULL);
   DBG_ASSERT(ma27b != NULL);
   DBG_ASSERT(ma27c != NULL);
   DBG_ASSERT(ma27i != NULL);

   options.GetNumericValue("ma27_pivtol", pivtol_, prefix);
   if( options.GetNumericValue("ma27_pivtolmax", pivtolmax_, prefix) )
   {
      ASSERT_EXCEPTION(pivtolmax_ >= pivtol_, OPTION_INVALID, "Option \"ma27_pivtolmax\": This value must be between "
                       "ma27_pivtol and 1.");
   }
   else
   {
      pivtolmax_ = Max(pivtolmax_, pivtol_);
   }

   Index print_level;
   options.GetIntegerValue("ma27_print_level", print_level, prefix);
   options.GetNumericValue("ma27_liw_init_factor", liw_init_factor_, prefix);
   options.GetNumericValue("ma27_la_init_factor", la_init_factor_, prefix);
   options.GetNumericValue("ma27_meminc_factor", meminc_factor_, prefix);
   options.GetBoolValue("ma27_skip_inertia_check", skip_inertia_check_, prefix);
   options.GetBoolValue("ma27_ignore_singularity", ignore_singularity_, prefix);
   // The following option is registered by OrigIpoptNLP
   options.GetBoolValue("warm_start_same_structure", warm_start_same_structure_, prefix);

   /* Set the default options for MA27 */
   ma27i(icntl_, cntl_);

   if( print_level == 0 )
   {
      icntl_[0] = 0;   // Suppress error messages
   }
   if( print_level <= 1 )
   {
      icntl_[1] = 0;   // Suppress warning messages
   }
   if( print_level >= 2 )
   {
      icntl_[2] = print_level - 2;   // diagnostic messages level
   }

   // Reset all private data
   initialized_ = false;
   pivtol_changed_ = false;
   refactorize_ = false;

   la_increase_ = false;
   liw_increase_ = false;

   if( !warm_start_same_structure_ )
   {
      dim_ = 0;
      nonzeros_ = 0;
   }
   else
   {
      ASSERT_EXCEPTION(dim_ > 0 && nonzeros_ > 0, INVALID_WARMSTART,
                       "Ma27TSolverInterface called with warm_start_same_structure, but the problem is solved for the first time.");
   }

   return true;
}

ESymSolverStatus Ma27TSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* airn,
   const Index* ajcn,
   Index        nrhs,
   Number*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("Ma27TSolverInterface::MultiSolve", dbg_verbosity);
   DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   DBG_ASSERT(initialized_);
   DBG_ASSERT(la_ != 0);

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

Number* Ma27TSolverInterface::GetValuesArrayPtr()
{
   DBG_START_METH("Ma27TSolverInterface::GetValuesArrayPtr", dbg_verbosity);
   DBG_ASSERT(initialized_);

   // If the size of a is to be increase for the next factorization
   // anyway, delete the current large array and just return enough
   // to store the values

   if( la_increase_ )
   {
      delete[] a_;
      a_ = NULL;
      a_ = new Number[nonzeros_];
   }

   return a_;
}

/** Initialize the local copy of the positions of the nonzero elements */
ESymSolverStatus Ma27TSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* airn,
   const Index* ajcn
)
{
   DBG_START_METH("Ma27TSolverInterface::InitializeStructure", dbg_verbosity);

   ESymSolverStatus retval = SYMSOLVER_SUCCESS;
   if( !warm_start_same_structure_ )
   {
      dim_ = dim;
      nonzeros_ = nonzeros;

      // Do the symbolic facotrization
      retval = SymbolicFactorization(airn, ajcn);
      if( retval != SYMSOLVER_SUCCESS )
      {
         return retval;
      }
   }
   else
   {
      ASSERT_EXCEPTION(dim_ == dim && nonzeros_ == nonzeros, INVALID_WARMSTART,
                       "Ma27TSolverInterface called with warm_start_same_structure, but the problem size has changed.");
   }

   initialized_ = true;

   return retval;
}

ESymSolverStatus Ma27TSolverInterface::SymbolicFactorization(
   const Index* airn,
   const Index* ajcn
)
{
   DBG_START_METH("Ma27TSolverInterface::SymbolicFactorization", dbg_verbosity);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
   }

   // Get memory for the IW workspace
   delete[] iw_;
   iw_ = NULL;

   // Overestimation factor for LIW (20% recommended in MA27 documentation)
   const Number LiwFact = 2.0;      // This is 100% overestimation
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "In Ma27TSolverInterface::InitializeStructure: Using overestimation factor LiwFact = %e\n", LiwFact);
   liw_ = (Index) (LiwFact * (Number(2 * nonzeros_ + 3 * dim_ + 1)));
   try
   {
      iw_ = new Index[liw_];
   }
   catch( const std::bad_alloc& )
   {
      Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Failed to allocate initial working space (iw_) for MA27\n");
      throw; // will be caught in IpIpoptApplication
   }

   // Get memory for IKEEP
   delete[] ikeep_;
   ikeep_ = NULL;
   ikeep_ = new Index[3 * dim_];

   if( Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA) )
   {
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                     "\nMatrix structure given to MA27 with dimension %" IPOPT_INDEX_FORMAT " and %" IPOPT_INDEX_FORMAT " nonzero entries:\n", dim_, nonzeros_);
      for( Index i = 0; i < nonzeros_; i++ )
      {
         Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA,
                        "A[%5d,%5d]\n", airn[i], ajcn[i]);
      }
   }

   // Call MA27AX
   Index N = dim_;
   Index NZ = nonzeros_;
   Index IFLAG = 0;
   Number OPS;
   Index INFO[20];
   Index* IW1 = new Index[2 * dim_];      // Get memory for IW1 (only local)
   ma27a(&N, &NZ, airn, ajcn, iw_, &liw_, ikeep_, IW1, &nsteps_, &IFLAG, icntl_, cntl_, INFO, &OPS);
   delete[] IW1;      // No longer required

   // Receive several information
   const Index& iflag = INFO[0];      // Information flag
   const Index& ierror = INFO[1];      // Error flag
   const Index& nrlnec = INFO[4];      // recommended value for la
   const Index& nirnec = INFO[5];      // recommended value for liw

   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Return values from MA27AD: IFLAG = %" IPOPT_INDEX_FORMAT ", IERROR = %" IPOPT_INDEX_FORMAT "\n", iflag, ierror);

   // Check if error occurred
   if( iflag != 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "*** Error from MA27AD *** IFLAG = %" IPOPT_INDEX_FORMAT " IERROR = %" IPOPT_INDEX_FORMAT "\n", iflag, ierror);
      if( iflag == 1 )
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "The index of a matrix is out of range.\nPlease check your implementation of the Jacobian and Hessian matrices.\n");
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
   }

   try
   {
      // Reserve memory for iw_ for later calls, based on suggested size
      delete[] iw_;
      iw_ = NULL;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Size of integer work space recommended by MA27 is %" IPOPT_INDEX_FORMAT "\n", nirnec);
      ComputeMemIncrease(liw_, liw_init_factor_ * (Number) nirnec, 0, "integer working space for MA27");
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Setting integer work space size to %" IPOPT_INDEX_FORMAT "\n", liw_);
      iw_ = new Index[liw_];

      // Reserve memory for a_
      delete[] a_;
      a_ = NULL;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Size of doublespace recommended by MA27 is %" IPOPT_INDEX_FORMAT "\n", nrlnec);
      ComputeMemIncrease(la_, la_init_factor_ * (Number) nrlnec, nonzeros_, "double working space for MA27");
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Setting double work space size to %" IPOPT_INDEX_FORMAT "\n", la_);
      a_ = new Number[la_];
   }
   catch( const std::bad_alloc& )
   {
      Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Failed to allocate more working space for MA27\n");
      throw; // will be caught in IpIpoptApplication
   }

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus Ma27TSolverInterface::Factorization(
   const Index* airn,
   const Index* ajcn,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("Ma27TSolverInterface::Factorization", dbg_verbosity);

   // Check if la should be increased
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().Start();
   }

   if( la_increase_ )
   {
      Number* a_old = a_;
      Index la_old = la_;
      ComputeMemIncrease(la_, meminc_factor_ * (Number) la_, 0, "double working space for MA27");
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: Increasing la from %" IPOPT_INDEX_FORMAT " to %" IPOPT_INDEX_FORMAT "\n", la_old, la_);
      try
      {
         a_ = new Number[la_];
      }
      catch( const std::bad_alloc& )
      {
         Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Failed to allocate more working space (a_) for MA27\n");
         throw; // will be caught in IpIpoptApplication
      }
      for( Index i = 0; i < nonzeros_; i++ )
      {
         a_[i] = a_old[i];
      }
      delete[] a_old;
      la_increase_ = false;
   }

   // Check if liw should be increased
   if( liw_increase_ )
   {
      delete[] iw_;
      iw_ = NULL;
      Index liw_old = liw_;
      ComputeMemIncrease(liw_, meminc_factor_ * (Number) liw_, 0, "integer working space for MA27");
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: Increasing liw from %" IPOPT_INDEX_FORMAT " to %" IPOPT_INDEX_FORMAT "\n", liw_old, liw_);
      try
      {
         iw_ = new Index[liw_];
      }
      catch( const std::bad_alloc& )
      {
         Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Failed to allocate more working space (iw_) for MA27\n");
         throw; // will be caught in IpIpoptApplication
      }
      liw_increase_ = false;
   }

   Index iflag;  // Information flag
   Index ncmpbr;  // Number of double precision compressions
   Index ncmpbi;  // Number of integer compressions

   // Call MA27BX; possibly repeatedly if workspaces are too small
   Index N = dim_;
   Index NZ = nonzeros_;
   Index* IW1 = new Index[2 * dim_];
   Index INFO[20];
   cntl_[0] = pivtol_;  // Set pivot tolerance

   ma27b(&N, &NZ, airn, ajcn, a_, &la_, iw_, &liw_, ikeep_, &nsteps_, &maxfrt_, IW1, icntl_, cntl_, INFO);
   delete[] IW1;

   // Receive information about the factorization
   iflag = INFO[0];  // Information flag
   const Index& ierror = INFO[1];  // Error flag
   ncmpbr = INFO[11];  // Number of double compressions
   ncmpbi = INFO[12];  // Number of integer compressions
   negevals_ = INFO[14];  // Number of negative eigenvalues

   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Return values from MA27BD: IFLAG = %" IPOPT_INDEX_FORMAT ", IERROR = %" IPOPT_INDEX_FORMAT "\n", iflag, ierror);

   DBG_PRINT((1, "Return from MA27BD iflag = %" IPOPT_INDEX_FORMAT " and ierror = %" IPOPT_INDEX_FORMAT "\n",
              iflag, ierror));

   // Check if factorization failed due to insufficient memory space
   // iflag==-3 if LIW too small (recommended value in ierror)
   // iflag==-4 if LA too small (recommended value in ierror)
   if( iflag == -3 || iflag == -4 )
   {
      // Increase size of both LIW and LA
      delete[] iw_;
      iw_ = NULL;
      delete[] a_;
      a_ = NULL;
      Index liw_old = liw_;
      Index la_old = la_;
      if( iflag == -3 )
      {
         ComputeMemIncrease(liw_, meminc_factor_ * (Number) ierror, 0, "integer working space for MA27");
         ComputeMemIncrease(la_, meminc_factor_ * (Number) la_, 0, "double working space for MA27");
      }
      else
      {
         ComputeMemIncrease(liw_, meminc_factor_ * (Number) liw_, 0, "integer working space for MA27");
         ComputeMemIncrease(la_, meminc_factor_ * (Number) ierror, 0, "double working space for MA27");
      }
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned iflag=%" IPOPT_INDEX_FORMAT " and requires more memory.\n Increase liw from %" IPOPT_INDEX_FORMAT " to %" IPOPT_INDEX_FORMAT " and la from %" IPOPT_INDEX_FORMAT " to %" IPOPT_INDEX_FORMAT " and factorize again.\n",
                     iflag, liw_old, liw_, la_old, la_);
      try
      {
         iw_ = new Index[liw_];
         a_ = new Number[la_];
      }
      catch( const std::bad_alloc& )
      {
         Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Failed to allocate more working space (iw_ and a_) for MA27\n");
         throw; // will be caught in IpIpoptApplication
      }
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_CALL_AGAIN;
   }

   // Check if the system is singular, and if some other error occurred
   if( iflag == -5 || (!ignore_singularity_ && iflag == 3) )
   {
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_SINGULAR;
   }
   else if( iflag == 3 )
   {
      Index missing_rank = dim_ - INFO[1];
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned iflag=%" IPOPT_INDEX_FORMAT " and detected rank deficiency of degree %" IPOPT_INDEX_FORMAT ".\n", iflag, missing_rank);
      // We correct the number of negative eigenvalues here to include
      // the zero eigenvalues, since otherwise we indicate the wrong
      // inertia.
      negevals_ += missing_rank;
   }
   else if( iflag != 0 )
   {
      // There is some error
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
   }

   // Check if it might be more efficient to use more memory next time
   // (if there were too many compressions for this factorization)
   if( ncmpbr >= 10 )
   {
      la_increase_ = true;
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned ncmpbr=%" IPOPT_INDEX_FORMAT ". Increase la before the next factorization.\n", ncmpbr);
   }
   if( ncmpbi >= 10 )
   {
      liw_increase_ = true;
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "MA27BD returned ncmpbi=%" IPOPT_INDEX_FORMAT ". Increase liw before the next factorization.\n", ncmpbr);
   }

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Number of doubles for MA27 to hold factorization (INFO(9)) = %" IPOPT_INDEX_FORMAT "\n", INFO[8]);
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Number of integers for MA27 to hold factorization (INFO(10)) = %" IPOPT_INDEX_FORMAT "\n", INFO[9]);

   // Check whether the number of negative eigenvalues matches the requested
   // count
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().End();
   }
   if( !skip_inertia_check_ && check_NegEVals && (numberOfNegEVals != negevals_) )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "In Ma27TSolverInterface::Factorization: negevals_ = %" IPOPT_INDEX_FORMAT ", but numberOfNegEVals = %" IPOPT_INDEX_FORMAT "\n", negevals_,
                     numberOfNegEVals);
      return SYMSOLVER_WRONG_INERTIA;
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus Ma27TSolverInterface::Backsolve(
   Index   nrhs,
   Number* rhs_vals
)
{
   DBG_START_METH("Ma27TSolverInterface::Backsolve", dbg_verbosity);
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().Start();
   }

   Index N = dim_;
   Number* W = new Number[maxfrt_];
   Index* IW1 = new Index[nsteps_];

   // For each right hand side, call MA27CX
   for( Index irhs = 0; irhs < nrhs; irhs++ )
   {
      if( DBG_VERBOSITY() >= 2 )
      {
         for( Index i = 0; i < dim_; i++ )
         {
            DBG_PRINT((2, "rhs[%5d] = %23.15e\n", i, rhs_vals[irhs * dim_ + i]));
         }
      }
      ma27c(&N, a_, &la_, iw_, &liw_, W, &maxfrt_, &rhs_vals[irhs * dim_], IW1, &nsteps_, icntl_, cntl_);
      if( DBG_VERBOSITY() >= 2 )
      {
         for( Index i = 0; i < dim_; i++ )
         {
            DBG_PRINT((2, "sol[%5d] = %23.15e\n", i, rhs_vals[irhs * dim_ + i]));
         }
      }
   }
   delete[] W;
   delete[] IW1;

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().End();
   }

   return SYMSOLVER_SUCCESS;
}

Index Ma27TSolverInterface::NumberOfNegEVals() const
{
   DBG_START_METH("Ma27TSolverInterface::NumberOfNegEVals", dbg_verbosity);
   DBG_ASSERT(ProvidesInertia());
   DBG_ASSERT(initialized_);
   return negevals_;
}

bool Ma27TSolverInterface::IncreaseQuality()
{
   DBG_START_METH("Ma27TSolverInterface::IncreaseQuality", dbg_verbosity);
   if( pivtol_ == pivtolmax_ )
   {
      return false;
   }

   pivtol_changed_ = true;

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Increasing pivot tolerance for MA27 from %7.2e ", pivtol_);
   pivtol_ = Min(pivtolmax_, std::pow(pivtol_, Number(0.75)));
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "to %7.2e.\n", pivtol_);
   return true;
}

} // namespace Ipopt
