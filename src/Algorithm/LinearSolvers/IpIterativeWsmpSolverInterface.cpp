// Copyright (C) 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter                 IBM    2009-09-18
//               based on IpWsmpSolverInterface.cpp (rev 1521)

#include "IpoptConfig.h"
#include "IpIterativeWsmpSolverInterface.hpp"
#include "IpBlas.hpp"

#include <cmath>
#include <cstdio>

/** Prototypes for WISMP's subroutines */
extern "C"
{
   void F77_FUNC(wsetmaxthrds, WSETMAXTHRDS)(
      const ipindex* NTHREADS
   );

   void F77_FUNC(wismp, WISMP)(
      const ipindex* N,
      const ipindex* IA,
      const ipindex* JA,
      const double*  AVALS,
      double*        B,
      const ipindex* LDB,
      double*        X,
      const ipindex* LDX,
      const ipindex* NRHS,
      double*        RMISC,
      double*        CVGH,
      ipindex*       IPARM,
      double*        DPARM
   );

   void F77_FUNC_(wsmp_clear, WSMP_CLEAR)(void);
}

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 3;
#endif

IterativeWsmpSolverInterface::IterativeWsmpSolverInterface()
   : a_(NULL),
     initialized_(false)
{
   DBG_START_METH("IterativeWsmpSolverInterface::IterativeWsmpSolverInterface()", dbg_verbosity);

   IPARM_ = new Index[64];
   DPARM_ = new double[64];
}

IterativeWsmpSolverInterface::~IterativeWsmpSolverInterface()
{
   DBG_START_METH("IterativeWsmpSolverInterface::~IterativeWsmpSolverInterface()",
                  dbg_verbosity);

   // Clear WSMP's memory
   F77_FUNC_(wsmp_clear, WSMP_CLEAR)();

   delete[] IPARM_;
   delete[] DPARM_;
   delete[] a_;
}

void IterativeWsmpSolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions)
{
   roptions->AddLowerBoundedIntegerOption(
      "wsmp_max_iter",
      "Maximal number of iterations in iterative WISMP", 1,
      1000,
      "",
      true);
   roptions->AddLowerBoundedNumberOption(
      "wsmp_inexact_droptol",
      "Drop tolerance for inexact factorization preconditioner in WISMP.",
      0.0, false,
      0.0,
      "DPARM(14) in WISMP",
      true);
   roptions->AddLowerBoundedNumberOption(
      "wsmp_inexact_fillin_limit",
      "Fill-in limit for inexact factorization preconditioner in WISMP.",
      0.0, false,
      0.0,
      "DPARM(15) in WISMP",
      true);
}

bool IterativeWsmpSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   options.GetIntegerValue("wsmp_num_threads", wsmp_num_threads_, prefix);

   Index wsmp_ordering_option;
   if( !options.GetIntegerValue("wsmp_ordering_option", wsmp_ordering_option, prefix) )
   {
      wsmp_ordering_option = 1;
   }

   Index wsmp_ordering_option2;
   if( !options.GetIntegerValue("wsmp_ordering_option2", wsmp_ordering_option2, prefix) )
   {
      wsmp_ordering_option2 = 0;
   }

   if( !options.GetNumericValue("wsmp_pivtol", wsmp_pivtol_, prefix) )
   {
      wsmp_pivtol_ = 1e-3;
   }

   if( options.GetNumericValue("wsmp_pivtolmax", wsmp_pivtolmax_, prefix) )
   {
      ASSERT_EXCEPTION(wsmp_pivtolmax_ >= wsmp_pivtol_, OPTION_INVALID,
                       "Option \"wsmp_pivtolmax\": This value must be between "
                       "wsmp_pivtol and 1.");
   }
   else
   {
      wsmp_pivtolmax_ = Max(wsmp_pivtolmax_, wsmp_pivtol_);
   }

   if( !options.GetIntegerValue("wsmp_scaling", wsmp_scaling_, prefix) )
   {
      wsmp_scaling_ = 1;
   }

   options.GetIntegerValue("wsmp_write_matrix_iteration", wsmp_write_matrix_iteration_, prefix);
   Index wsmp_max_iter;
   options.GetIntegerValue("wsmp_max_iter", wsmp_max_iter, prefix);
   options.GetNumericValue("wsmp_inexact_droptol", wsmp_inexact_droptol_, prefix);
   options.GetNumericValue("wsmp_inexact_fillin_limit", wsmp_inexact_fillin_limit_, prefix);

   // Reset all private data
   dim_ = 0;
   initialized_ = false;
   pivtol_changed_ = false;
   have_symbolic_factorization_ = false;
   delete[] a_;
   a_ = NULL;

#if 1
   // Set the number of threads
   Index NTHREADS = wsmp_num_threads_;
   F77_FUNC(wsetmaxthrds, WSETMAXTHRDS)(&NTHREADS);
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "WSMP will use %" IPOPT_INDEX_FORMAT " threads.\n", wsmp_num_threads_);
#else
   Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                  "Not setting WISMP threads at the moment.\n");
#endif

   // Get WSMP's default parameters and set the ones we want differently
   IPARM_[0] = 0;
   IPARM_[1] = 0;
   IPARM_[2] = 0;
   Index idmy = 0;
   double ddmy = 0.;
   F77_FUNC(wismp, WISMP)(&idmy, &idmy, &idmy, &ddmy, &ddmy, &idmy, &ddmy, &idmy, &idmy, &ddmy, &ddmy, IPARM_, DPARM_);

   if( IPARM_[63] < 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA, "Error %" IPOPT_INDEX_FORMAT " from WSMP initialization.\n", IPARM_[63]);
      return false;
   }

   IPARM_[3] = 3; // Upper triangular portion of matrix in CSR format
   // (same as for WSSMP)
   IPARM_[6] = 3;
   IPARM_[13] = 0; // do not overwrite avals
   IPARM_[27] = 0; // to make runs repeatable

#if 1
   IPARM_[5] = wsmp_max_iter; // maximal number of iterations
   IPARM_[15] = wsmp_ordering_option; // ordering option
   IPARM_[16] = wsmp_ordering_option2; // for ordering in IP methods?
#endif
   DPARM_[13] = wsmp_inexact_droptol_;
   DPARM_[14] = wsmp_inexact_fillin_limit_;

   // DELETE
   IPARM_[33] = 0;

   matrix_file_number_ = 0;

   // Check for SPINLOOPTIME and YIELDLOOPTIME?

   return true;
}

ESymSolverStatus IterativeWsmpSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   double*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("IterativeWsmpSolverInterface::MultiSolve", dbg_verbosity);
   DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   DBG_ASSERT(initialized_);

   // check if a factorization has to be done
   if( new_matrix || pivtol_changed_ )
   {
      pivtol_changed_ = false;
      // perform the factorization
      ESymSolverStatus retval;
      retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
      if( retval != SYMSOLVER_SUCCESS )
      {
         DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
         return retval;  // Matrix singular or error occurred
      }
   }

   // do the solve
   return Solve(ia, ja, nrhs, rhs_vals);
}

double* IterativeWsmpSolverInterface::GetValuesArrayPtr()
{
   DBG_ASSERT(initialized_);
   DBG_ASSERT(a_);
   return a_;
}

/* Initialize the local copy of the positions of the nonzero elements */
ESymSolverStatus IterativeWsmpSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* ia,
   const Index* ja
)
{
   DBG_START_METH("IterativeWsmpSolverInterface::InitializeStructure", dbg_verbosity);
   dim_ = dim;

   // Make space for storing the matrix elements
   delete[] a_;
   a_ = NULL;
   a_ = new double[nonzeros];

   // Do the symbolic factorization
   ESymSolverStatus retval = SymbolicFactorization(ia, ja);
   if( retval != SYMSOLVER_SUCCESS )
   {
      return retval;
   }

   initialized_ = true;

   return retval;
}

ESymSolverStatus IterativeWsmpSolverInterface::SymbolicFactorization(
   const Index* /*ia*/,
   const Index* /*ja*/
)
{
   DBG_START_METH("IterativeWsmpSolverInterface::SymbolicFactorization",
                  dbg_verbosity);

   // This is postponed until the first factorization call, since
   // then the values in the matrix are known
   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus IterativeWsmpSolverInterface::InternalSymFact(
   const Index* ia,
   const Index* ja
)
{
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
   }

   // Call WISMP for ordering and symbolic factorization
   Index N = dim_;
   IPARM_[1] = 1; // ordering
   IPARM_[2] = 1; // symbolic factorization
   Index idmy = 0;
   double ddmy = 0.;
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Calling WISMP-1-1 for symbolic analysis.\n");
   F77_FUNC(wismp, WISMP)(&N, ia, ja, a_, &ddmy, &idmy, &ddmy, &idmy, &idmy, &ddmy, &ddmy, IPARM_, DPARM_);
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Done with WISMP-1-1 for symbolic analysis.\n");

   Index ierror = IPARM_[63];
   if( ierror != 0 )
   {
      if( ierror == -102 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error: WISMP is not able to allocate sufficient amount of memory during ordering/symbolic factorization.\n");
      }
      else if( ierror > 0 )
      {
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Matrix appears to be singular (with ierror = %" IPOPT_INDEX_FORMAT ").\n", ierror);
         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemSymbolicFactorization().End();
         }
         return SYMSOLVER_SINGULAR;
      }
      else
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error in WISMP during ordering/symbolic factorization phase.\n     Error code is %" IPOPT_INDEX_FORMAT ".\n", ierror);
      }
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
   }
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Predicted memory usage for WISMP after symbolic factorization IPARM(23)= %" IPOPT_INDEX_FORMAT ".\n", IPARM_[22]);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus IterativeWsmpSolverInterface::Factorization(
   const Index* ia,
   const Index* ja,
   bool         /*check_NegEVals*/,
   Index        /*numberOfNegEVals*/
)
{
   DBG_START_METH("IterativeWsmpSolverInterface::Factorization", dbg_verbosity);

   // If desired, write out the matrix
   Index iter_count = -1;
   if( HaveIpData() )
   {
      iter_count = IpData().iter_count();
   }
   if( iter_count == wsmp_write_matrix_iteration_ )
   {
      matrix_file_number_++;
      char buf[256];
      Snprintf(buf, 255, "wsmp_matrix_%" IPOPT_INDEX_FORMAT "_%" IPOPT_INDEX_FORMAT ".dat", iter_count, matrix_file_number_);
      Jnlst().Printf(J_SUMMARY, J_LINEAR_ALGEBRA,
                     "Writing WSMP matrix into file %s.\n", buf);
      FILE* fp = fopen(buf, "w");
      fprintf(fp, "%" IPOPT_INDEX_FORMAT "\n", dim_); // N
      for( Index icol = 0; icol < dim_; icol++ )
      {
         fprintf(fp, "%" IPOPT_INDEX_FORMAT "", ia[icol + 1] - ia[icol]); // number of elements for this column
         // Now for each column we write row indices and values
         for( Index irow = ia[icol]; irow < ia[icol + 1]; irow++ )
         {
            fprintf(fp, " %23.16e %" IPOPT_INDEX_FORMAT "", a_[irow - 1], ja[irow - 1]);
         }
         fprintf(fp, "\n");
      }
      fclose(fp);
   }

   // Check if we have to do the symbolic factorization and ordering
   // phase yet
   if( !have_symbolic_factorization_ )
   {
      ESymSolverStatus retval = InternalSymFact(ia, ja);
      if( retval != SYMSOLVER_SUCCESS )
      {
         return retval;
      }
      have_symbolic_factorization_ = true;
   }

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().Start();
   }

   // Call WSSMP for numerical factorization
   Index N = dim_;
   IPARM_[1] = 2; // value analysis
   IPARM_[2] = 3; // preconditioner generation
   DPARM_[10] = wsmp_pivtol_; // set current pivot tolerance
   Index idmy = 0;
   double ddmy = 0.;

   // set drop tolerances for now....
   if( wsmp_inexact_droptol_ != 0. )
   {
      DPARM_[13] = wsmp_inexact_droptol_;
   }
   if( wsmp_inexact_fillin_limit_ != 0. )
   {
      DPARM_[14] = wsmp_inexact_fillin_limit_;
   }

   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Calling WISMP-2-3 with DPARM(14) = %8.2e and DPARM(15) = %8.2e.\n", DPARM_[13], DPARM_[14]);
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Calling WISMP-2-3 for value analysis and preconditioner computation.\n");
   F77_FUNC(wismp, WISMP)(&N, ia, ja, a_, &ddmy, &idmy, &ddmy, &idmy, &idmy, &ddmy, &ddmy, IPARM_, DPARM_);
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Done with WISMP-2-3 for value analysis and preconditioner computation.\n");
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Done with WISMP-2-3 with DPARM(14) = %8.2e and DPARM(15) = %8.2e.\n", DPARM_[13], DPARM_[14]);
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "                         DPARM(4) = %8.2e and DPARM(5) = %8.2e and ratio = %8.2e.\n", DPARM_[3], DPARM_[4],
                  DPARM_[3] / DPARM_[4]);

   const Index ierror = IPARM_[63];
   if( ierror > 0 )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "WISMP detected that the matrix is singular and encountered %" IPOPT_INDEX_FORMAT " zero pivots.\n", dim_ + 1 - ierror);
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_SINGULAR;
   }
   else if( ierror != 0 )
   {
      if( ierror == -102 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error: WISMP is not able to allocate sufficient amount of memory during factorization.\n");
      }
      else
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error in WSMP during factorization phase.\n     Error code is %" IPOPT_INDEX_FORMAT ".\n", ierror);
      }
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_FATAL_ERROR;
   }
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Memory usage for WISMP after factorization IPARM(23) = %" IPOPT_INDEX_FORMAT "\n", IPARM_[22]);

#if 0
   // Check whether the number of negative eigenvalues matches the requested
   // count
   if (check_NegEVals && (numberOfNegEVals != negevals_))
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %" IPOPT_INDEX_FORMAT ", but we got %" IPOPT_INDEX_FORMAT ".\n",
                     numberOfNegEVals, negevals_);
      if (HaveIpData())
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      return SYMSOLVER_WRONG_INERTIA;
   }
#endif

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemFactorization().End();
   }
   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus IterativeWsmpSolverInterface::Solve(
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   double*      rhs_vals
)
{
   DBG_START_METH("IterativeWsmpSolverInterface::Solve", dbg_verbosity);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().Start();
   }

   // Call WISMP to solve for some right hand sides.  The solution
   // will be stored in rhs_vals, and we need to make a copy of the
   // original right hand side before the call.
   Index N = dim_;
   Index LDB = dim_;
   double* RHS = new double[dim_ * nrhs];
   IpBlasCopy(dim_ * nrhs, rhs_vals, 1, RHS, 1);
   Index LDX = dim_; // Q: Do we have to zero out solution?
   Index NRHS = nrhs;
   IPARM_[1] = 4; // Iterative solver solution
   IPARM_[2] = 4;

   double* CVGH = NULL;
   if( Jnlst().ProduceOutput(J_MOREDETAILED, J_LINEAR_ALGEBRA) )
   {
      IPARM_[26] = 1; // Record convergence history
      CVGH = new double[IPARM_[5] + 1];
   }
   else
   {
      IPARM_[26] = 0;
   }

   double ddmy;
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Calling WISMP-4-4 for backsolve.\n");
   F77_FUNC(wismp, WISMP)(&N, ia, ja, a_, RHS, &LDB, rhs_vals, &LDX, &NRHS, &ddmy, CVGH, IPARM_, DPARM_);
   Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                  "Done with WISMP-4-4 for backsolve.\n");
   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().End();
   }

   Index ierror = IPARM_[63];
   if( ierror != 0 )
   {
      if( ierror == -102 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error: WISMP is not able to allocate sufficient amount of memory during ordering/symbolic factorization.\n");
      }
      else
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error in WISMP during ordering/symbolic factorization phase.\n     Error code is %" IPOPT_INDEX_FORMAT ".\n", ierror);
      }
      delete[] CVGH;
      return SYMSOLVER_FATAL_ERROR;
   }
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Number of iterative solver steps in WISMP: %" IPOPT_INDEX_FORMAT "\n", IPARM_[25]);
   if( Jnlst().ProduceOutput(J_MOREDETAILED, J_LINEAR_ALGEBRA) )
   {
      DBG_ASSERT(CVGH != NULL);
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "WISMP congergence history:\n");
      for( Index i = 0; i <= IPARM_[25]; ++i )
      {
         Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                        " Resid[%3d] = %13.6e\n", i, CVGH[i]);
      }
   }
   delete[] CVGH;

   return SYMSOLVER_SUCCESS;
}

Index IterativeWsmpSolverInterface::NumberOfNegEVals() const
{
   DBG_START_METH("IterativeWsmpSolverInterface::NumberOfNegEVals", dbg_verbosity);
   DBG_ASSERT(false);
   return -1;
}

bool IterativeWsmpSolverInterface::IncreaseQuality()
{
   // TODO!!!!
   return false;
#if 0
   DBG_START_METH("IterativeWsmpSolverInterface::IncreaseQuality", dbg_verbosity);
   if( wsmp_pivtol_ == wsmp_pivtolmax_ )
   {
      return false;
   }
   pivtol_changed_ = true;

   if( wsmp_inexact_droptol_ != 0. )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Increasing dropotol for WISMP from %7.2e ", wsmp_inexact_droptol_);
      wsmp_inexact_droptol_ = DPARM_[13];
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "to %7.2e (suggested value).\n", wsmp_inexact_droptol_);
   }
   else
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Not increasing dropotol for WISMP, it is just reusing new value");
   }
   if( wsmp_inexact_fillin_limit_ != 0. )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Increasing fillin limit for WISMP from %7.2e ", wsmp_inexact_fillin_limit_);
      wsmp_inexact_fillin_limit_ = DPARM_[14];
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "to %7.2e (suggested value).\n", wsmp_inexact_fillin_limit_);
   }
   else
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Not increasing fillin limit for WISMP, it is just reusing new value");
   }

   return true;
#endif
}

} // namespace Ipopt
