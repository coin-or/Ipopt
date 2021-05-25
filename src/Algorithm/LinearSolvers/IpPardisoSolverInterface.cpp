// Copyright (C) 2005, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17
//
//           Olaf Schenk                      Univ of Basel 2005-09-20
//                  - changed options, added PHASE_ flag

#include "IpoptConfig.h"
#include "IpPardisoSolverInterface.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <fstream>
#include <iomanip>

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

/** pointer to Pardiso function that can be set via PardisoSolverInterface::SetFunctions() */
static IPOPT_DECL_PARDISOINIT(*user_pardisoinit) = NULL;
static IPOPT_DECL_PARDISO(*user_pardiso) = NULL;
static bool user_isparallel = false;
#ifdef PARDISO_MATCHING_PREPROCESS
static IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*user_smat_reordering_pardiso_wsmp) = NULL;
#endif

PardisoSolverInterface::PardisoSolverInterface(
   SmartPtr<LibraryLoader> pardisoloader_
)  : a_(NULL),
#ifdef PARDISO_MATCHING_PREPROCESS
   ia2(NULL),
   ja2(NULL),
   a2_(NULL),
   perm2(NULL),
   scale2(NULL),
#endif
   negevals_(-1),
   initialized_(false),
   MAXFCT_(1),
   MNUM_(1),
   MTYPE_(-2),
   MSGLVL_(0),
   debug_last_iter_(-1),
   pardisoloader(pardisoloader_),
   pardisoinit(NULL),
   pardiso(NULL),
#ifdef PARDISO_MATCHING_PREPROCESS
   smat_reordering_pardiso_wsmp(NULL),
#endif
   pardiso_exist_parallel(false)
{
   DBG_START_METH("PardisoSolverInterface::PardisoSolverInterface()", dbg_verbosity);

   PT_ = new void* [64];
   IPARM_ = new Index[64];
   DPARM_ = new Number[64];
}

PardisoSolverInterface::~PardisoSolverInterface()
{
   DBG_START_METH("PardisoSolverInterface::~PardisoSolverInterface()",
                  dbg_verbosity);

   // Tell Pardiso to release all memory
   if( initialized_ )
   {
      Index PHASE = -1;
      Index N = dim_;
      Index NRHS = 0;
      Index ERROR;
      Index idmy = 0;
      Number ddmy = 0.;
      pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N, &ddmy, &idmy, &idmy, &idmy, &NRHS, IPARM_, &MSGLVL_, &ddmy,
              &ddmy, &ERROR, DPARM_);
      DBG_ASSERT(ERROR == 0);
   }

   delete[] PT_;
   delete[] IPARM_;
   delete[] DPARM_;
   delete[] a_;

#ifdef PARDISO_MATCHING_PREPROCESS
   delete[] ia2;
   delete[] ja2;
   delete[] a2_;
   delete[] perm2;
   delete[] scale2;
#endif

}

void PardisoSolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddStringOption3(
      "pardiso_matching_strategy",
      "Matching strategy to be used by Pardiso",
      "complete+2x2",
      "complete", "Match complete (IPAR(13)=1)",
      "complete+2x2", "Match complete+2x2 (IPAR(13)=2)",
      "constraints", "Match constraints (IPAR(13)=3)",
      "This is IPAR(13) in Pardiso manual.");
   roptions->AddStringOption2(
      "pardiso_redo_symbolic_fact_only_if_inertia_wrong",
      "Toggle for handling case when elements were perturbed by Pardiso.",
      "no",
      "no", "Always redo symbolic factorization when elements were perturbed",
      "yes", "Only redo symbolic factorization when elements were perturbed if also the inertia was wrong",
      "",
      true);
   roptions->AddBoolOption(
      "pardiso_repeated_perturbation_means_singular",
      "Whether to assume that matrix is singular if elements were perturbed after recent symbolic factorization.",
      false,
      "",
      true);
   //roptions->AddLowerBoundedIntegerOption(
   //  "pardiso_out_of_core_power",
   //  "Enables out-of-core variant of Pardiso",
   //  0, 0,
   //  "Setting this option to a positive integer k makes Pardiso work in the "
   //  "out-of-core variant where the factor is split in 2^k subdomains.  This "
   //  "is IPARM(50) in the Pardiso manual.  This option is only available if "
   //  "Ipopt has been compiled with Pardiso.");
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_msglvl",
      "Pardiso message level",
      0,
      0,
      "This is MSGLVL in the Pardiso manual.");
   roptions->AddBoolOption(
      "pardiso_skip_inertia_check",
      "Whether to pretend that inertia is correct.",
      false,
      "Setting this option to \"yes\" essentially disables inertia check. "
      "This option makes the algorithm non-robust and easily fail, but it might give some insight into the necessity of inertia control.",
      true);
   roptions->AddIntegerOption(
      "pardiso_max_iterative_refinement_steps",
      "Limit on number of iterative refinement steps.",
      0,
      "The solver does not perform more than the absolute value of this value steps of iterative refinement and "
      "stops the process if a satisfactory level of accuracy of the solution in terms of backward error is achieved. "
      "If negative, the accumulation of the residue uses extended precision real and complex data types. "
      "Perturbed pivots result in iterative refinement. "
      "The solver automatically performs two steps of iterative refinements when perturbed pivots are obtained during the numerical factorization and this option is set to 0.");
   roptions->AddStringOption6(
      "pardiso_order",
      "Controls the fill-in reduction ordering algorithm for the input matrix.",
      "metis",
      "amd", "minimum degree algorithm",
      "one", "",
      "metis", "MeTiS nested dissection algorithm",
      "pmetis", "parallel (OpenMP) version of MeTiS nested dissection algorithm",
      "four", "",
      "five", "");
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_max_iter",
      "Maximum number of Krylov-Subspace Iteration",
      1,
      500,
      "DPARM(1)",
      true);
   roptions->AddBoundedNumberOption(
      "pardiso_iter_relative_tol",
      "Relative Residual Convergence",
      0.0, true,
      1.0, true,
      1e-6,
      "DPARM(2)",
      true);
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_iter_coarse_size",
      "Maximum Size of Coarse Grid Matrix",
      1,
      5000,
      "DPARM(3)",
      true);
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_iter_max_levels",
      "Maximum Size of Grid Levels",
      1,
      10,
      "DPARM(4)",
      true);
   roptions->AddBoundedNumberOption(
      "pardiso_iter_dropping_factor",
      "dropping value for incomplete factor",
      0.0, true,
      1.0, true,
      0.5,
      "DPARM(5)",
      true);
   roptions->AddBoundedNumberOption(
      "pardiso_iter_dropping_schur",
      "dropping value for sparsify schur complement factor",
      0.0, true,
      1.0, true,
      1e-1,
      "DPARM(6)",
      true);
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_iter_max_row_fill",
      "max fill for each row",
      1,
      10000000,
      "DPARM(7)",
      true);
   roptions->AddLowerBoundedNumberOption(
      "pardiso_iter_inverse_norm_factor",
      "",
      1, true,
      5000000,
      "DPARM(8)",
      true);
   roptions->AddBoolOption(
      "pardiso_iterative",
      "Switch for iterative solver in Pardiso library",
      false,
      "",
      true);
   roptions->AddLowerBoundedIntegerOption(
      "pardiso_max_droptol_corrections",
      "Maximal number of decreases of drop tolerance during one solve.",
      1,
      4,
      "This is relevant only for iterative Pardiso options.",
      true);
}

void PardisoSolverInterface::SetFunctions(
   IPOPT_DECL_PARDISOINIT(*pardisoinit),
   IPOPT_DECL_PARDISO(*pardiso),
   bool isparallel,
   IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*smat_reordering_pardiso_wsmp)
)
{
   DBG_ASSERT(pardisoinit != NULL);
   DBG_ASSERT(pardiso != NULL);
#ifdef PARDISO_MATCHING_PREPROCESS
   DBG_ASSERT(smat_reordering_pardiso_wsmp != NULL);
#endif

   user_pardisoinit = pardisoinit;
   user_pardiso = pardiso;
   user_isparallel = isparallel;
#ifdef PARDISO_MATCHING_PREPROCESS
   user_smat_reordering_pardiso_wsmp = smat_reordering_pardiso_wsmp;
#else
   (void)smat_reordering_pardiso_wsmp;
#endif
}

bool PardisoSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( user_pardisoinit != NULL )
   {
      pardisoinit = user_pardisoinit;
      pardiso = user_pardiso;
      pardiso_exist_parallel = user_isparallel;
#ifdef PARDISO_MATCHING_PREPROCESS
      smat_reordering_pardiso_wsmp = user_smat_reordering_pardiso_wsmp;
#endif
   }
   else
   {
      DBG_ASSERT(IsValid(pardisoloader));

      pardisoinit = (IPOPT_DECL_PARDISOINIT(*))pardisoloader->loadSymbol("pardisoinit");
      pardiso = (IPOPT_DECL_PARDISO(*))pardisoloader->loadSymbol("pardiso");
#ifdef PARDISO_MATCHING_PREPROCESS
      smat_reordering_pardiso_wsmp = (IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*))pardisoloader->loadSymbol("smat_reordering_pardiso_wsmp");
#endif
      // load pardiso_ipopt_newinterface only as check that we get Pardiso >= 4.0.0 from pardiso-project
      pardisoloader->loadSymbol("pardiso_ipopt_newinterface");
      // set pardiso_exist_parallel to true if symbol pardiso_exist_parallel exists in pardiso library
      // ignore if symbol doesn't exist (pardiso_exist_parallel == false)
      try
      {
         pardisoloader->loadSymbol("pardiso_exist_parallel");
         pardiso_exist_parallel = true;
      }
      catch( const DYNAMIC_LIBRARY_FAILURE& )
      {
         DBG_ASSERT(!pardiso_exist_parallel);
      }
   }

   DBG_ASSERT(pardisoinit != NULL);
   DBG_ASSERT(pardiso != NULL);
#ifdef PARDISO_MATCHING_PREPROCESS
   DBG_ASSERT(smat_reordering_pardiso_wsmp);
#endif

   Index enum_int;
   options.GetEnumValue("pardiso_matching_strategy", enum_int, prefix);
   match_strat_ = PardisoMatchingStrategy(enum_int);
   options.GetBoolValue("pardiso_redo_symbolic_fact_only_if_inertia_wrong",
                        pardiso_redo_symbolic_fact_only_if_inertia_wrong_, prefix);
   options.GetBoolValue("pardiso_repeated_perturbation_means_singular", pardiso_repeated_perturbation_means_singular_,
                        prefix);
   //Index pardiso_out_of_core_power;
   //options.GetIntegerValue("pardiso_out_of_core_power",
   //                        pardiso_out_of_core_power, prefix);
   options.GetBoolValue("pardiso_skip_inertia_check", skip_inertia_check_, prefix);
   Index pardiso_msglvl;
   options.GetIntegerValue("pardiso_msglvl", pardiso_msglvl, prefix);
   Index max_iterref_steps;
   options.GetIntegerValue("pardiso_max_iterative_refinement_steps", max_iterref_steps, prefix);
   Index order;
   options.GetEnumValue("pardiso_order", order, prefix);
   options.GetBoolValue("pardiso_iterative", pardiso_iterative_, prefix);
   Index pardiso_max_iter;
   options.GetIntegerValue("pardiso_max_iter", pardiso_max_iter, prefix);
   Number pardiso_iter_relative_tol;
   options.GetNumericValue("pardiso_iter_relative_tol", pardiso_iter_relative_tol, prefix);
   Index pardiso_iter_coarse_size;
   options.GetIntegerValue("pardiso_iter_coarse_size", pardiso_iter_coarse_size, prefix);
   Index pardiso_iter_max_levels;
   options.GetIntegerValue("pardiso_iter_max_levels", pardiso_iter_max_levels, prefix);
   Number pardiso_iter_dropping_factor;
   options.GetNumericValue("pardiso_iter_dropping_factor", pardiso_iter_dropping_factor, prefix);
   Number pardiso_iter_dropping_schur;
   options.GetNumericValue("pardiso_iter_dropping_schur", pardiso_iter_dropping_schur, prefix);
   Index pardiso_iter_max_row_fill;
   options.GetIntegerValue("pardiso_iter_max_row_fill", pardiso_iter_max_row_fill, prefix);
   Number pardiso_iter_inverse_norm_factor;
   options.GetNumericValue("pardiso_iter_inverse_norm_factor", pardiso_iter_inverse_norm_factor, prefix);
   options.GetIntegerValue("pardiso_max_droptol_corrections", pardiso_max_droptol_corrections_, prefix);

   // Number value = 0.0;

   // Tell Pardiso to release all memory if it had been used before
   if( initialized_ )
   {
      Index PHASE = -1;
      Index N = dim_;
      Index NRHS = 0;
      Index ERROR;
      Index idmy;
      Number ddmy;
      pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N, &ddmy, &idmy, &idmy, &idmy, &NRHS, IPARM_, &MSGLVL_, &ddmy,
              &ddmy, &ERROR, DPARM_);
      DBG_ASSERT(ERROR == 0);
   }

   // Reset all private data
   dim_ = 0;
   nonzeros_ = 0;
   have_symbolic_factorization_ = false;
   initialized_ = false;
   delete[] a_;
   a_ = NULL;

#ifdef PARDISO_MATCHING_PREPROCESS
   delete[] ia2;
   ia2 = NULL;

   delete[] ja2;
   ja2 = NULL;

   delete[] a2_;
   a2_ = NULL;

   delete[] perm2;
   perm2 = NULL;

   delete[] scale2;
   scale2 = NULL;
#endif

   // Call Pardiso's initialization routine
   memset(PT_, 0, 64); // needs to be initialized to 0 according to MKL Pardiso docu
   IPARM_[0] = 0;  // Tell it to fill IPARM with default values(?)

   Index ERROR = 0;
   Index SOLVER = 0; // initialize only direct solver

   pardisoinit(PT_, &MTYPE_, &SOLVER, IPARM_, DPARM_, &ERROR);

   if( ERROR != 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA, "Problem with Pardiso license (error %" IPOPT_INDEX_FORMAT " from pardisoinit).\n", ERROR);
      return false;
   }

   // Set some parameters for Pardiso
   IPARM_[0] = 1;  // Don't use the default values

   int num_procs = 1;
   if( pardiso_exist_parallel )
   {
      // Obtain the numbers of processors from the value of OMP_NUM_THREADS
      char* var = getenv("OMP_NUM_THREADS");
      if( var != NULL )
      {
         sscanf(var, "%d", &num_procs);
         if( num_procs < 1 )
         {
            Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                           "Invalid value for OMP_NUM_THREADS (\"%s\").\n", var);
            return false;
         }
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Using environment OMP_NUM_THREADS = %d as the number of processors for PARDISO.\n", num_procs);
      }
   }
   else
   {
      Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                     "You should set the environment variable OMP_NUM_THREADS to the number of processors used in Pardiso (e.g., 1).\n\n");
   }

   IPARM_[1] = order;
   IPARM_[2] = num_procs; // Set the number of processors
   IPARM_[5] = 1;  // Overwrite right-hand side
   IPARM_[7] = max_iterref_steps;

   // Options suggested by Olaf Schenk
   IPARM_[9] = 12;
   IPARM_[10] = 2; // Results in better scaling
   // Matching information:  IPARM_[12] = 1 seems ok, but results in a
   // large number of pivot perturbation
   // Matching information:  IPARM_[12] = 2 robust,  but more  expensive method
   IPARM_[12] = (int) match_strat_;

   IPARM_[20] = 3; // Results in better accuracy
   IPARM_[23] = 1; // parallel fac
   IPARM_[24] = 1; // parallel solve
#ifdef IPOPT_SINGLE
   IPARM_[28] = 1; // 32-bit factorization (single precision)
#else
   IPARM_[28] = 0; // 64-bit factorization (double precision)
#endif
   IPARM_[29] = 80; // we need this for IPOPT interface
   //IPARM_[33] = 1; // bit-by-bit identical results in parallel run

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Pardiso matrix ordering     (IPARM(2)): %" IPOPT_INDEX_FORMAT "\n", IPARM_[1]);
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Pardiso max. iterref. steps (IPARM(8)): %" IPOPT_INDEX_FORMAT "\n", IPARM_[7]);
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Pardiso matching strategy  (IPARM(13)): %" IPOPT_INDEX_FORMAT "\n", IPARM_[12]);

   if( pardiso_iterative_ )
   {
      IPARM_[31] = 1;  // active direct solver

      DPARM_[0] = pardiso_max_iter; // maximum number of Krylov-Subspace Iteration
      // Default is 300
      // 1 <= value <= e.g. 1000
      DPARM_[1] = pardiso_iter_relative_tol; // Relative Residual Convergence
      // e.g.  pardiso_iter_tol
      // Default is 1e-6
      // 1e-16 <= value < 1
      DPARM_[2] = pardiso_iter_coarse_size; // Maximum Size of Coarse Grid Matrix
      // e.g.  pardiso_coarse_grid
      // Default is 5000
      // 1 <= value < number of equations
      DPARM_[3] = pardiso_iter_max_levels; // Maximum Number of Grid Levels
      // e.g.  pardiso_max_grid
      // Default is 10000
      // 1 <= value < number of equations
      DPARM_[4] = pardiso_iter_dropping_factor;  // dropping value for incomplete factor
      // e.g.  pardiso_dropping_factor
      // Default is 0.5
      // 1e-16 <= value < 1
      DPARM_[5] = pardiso_iter_dropping_schur;  // dropping value for sparsify schur complementfactor
      // e.g.  pardiso_dropping_schur
      // Default is 0.1
      // 1e-16 <= value < 1
      DPARM_[6] = pardiso_iter_max_row_fill;  // max fill for each row
      // e.g.  pardiso_max_fill
      // Default is 1000
      // 1 <= value < 100000
      DPARM_[7] = pardiso_iter_inverse_norm_factor;  // dropping value for sparsify schur complementfactor
      // e.g.  pardiso_inverse_norm_factor
      // Default is 500
      // 2 <= value < 50000
      DPARM_[8] = 25; // maximum number of non-improvement steps
   }

   MSGLVL_ = pardiso_msglvl;

   // Option for the out of core variant
   //IPARM_[49] = pardiso_out_of_core_power;

   return true;
}

ESymSolverStatus PardisoSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   Number*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("PardisoSolverInterface::MultiSolve", dbg_verbosity);
   DBG_ASSERT(!check_NegEVals || ProvidesInertia());
   DBG_ASSERT(initialized_);

   // check if a factorization has to be done
   if( new_matrix )
   {
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

Number* PardisoSolverInterface::GetValuesArrayPtr()
{
   DBG_ASSERT(initialized_);
   DBG_ASSERT(a_);
   return a_;
}

ESymSolverStatus PardisoSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* ia,
   const Index* ja
)
{
   DBG_START_METH("PardisoSolverInterface::InitializeStructure", dbg_verbosity);
   dim_ = dim;
   nonzeros_ = nonzeros;

   // Make space for storing the matrix elements
   delete[] a_;
   a_ = NULL;
   a_ = new Number[nonzeros_];

   // Do the symbolic facotrization
   ESymSolverStatus retval = SymbolicFactorization(ia, ja);
   if( retval != SYMSOLVER_SUCCESS )
   {
      return retval;
   }

   initialized_ = true;

   return retval;
}

ESymSolverStatus PardisoSolverInterface::SymbolicFactorization(
   const Index* /*ia*/,
   const Index* /*ja*/
)
{
   DBG_START_METH("PardisoSolverInterface::SymbolicFactorization",
                  dbg_verbosity);

   // Since Pardiso requires the values of the nonzeros of the matrix
   // for an efficient symbolic factorization, we postpone that task
   // until the first call of Factorize.  All we do here is to reset
   // the flag (in case this interface is called for a matrix with a
   // new structure).

   have_symbolic_factorization_ = false;

   return SYMSOLVER_SUCCESS;
}

static
void write_iajaa_matrix(
   Index        N,
   const Index* ia,
   const Index* ja,
   Number*      a_,
   Number*      rhs_vals,
   int          iter_cnt,
   int          sol_cnt
)
{
   if( getenv("IPOPT_WRITE_MAT") )
   {
      /* Write header */
      char mat_name[128];
      char mat_pref[32];

      Index NNZ = ia[N] - 1;
      Index i;

      if( getenv("IPOPT_WRITE_PREFIX") )
      {
         strcpy(mat_pref, getenv("IPOPT_WRITE_PREFIX"));
      }
      else
      {
         strcpy(mat_pref, "mat-ipopt");
      }

      Snprintf(mat_name, 127, "%s_%03d-%02d.iajaa", mat_pref, iter_cnt, sol_cnt);

      // Open and write matrix file.
      std::ofstream mat_file(mat_name);
      mat_file << std::setprecision(std::numeric_limits<Number>::digits10 + 1);

      mat_file << N << std::endl;
      mat_file << NNZ << std::endl;

      for( i = 0; i < N + 1; i++ )
      {
         mat_file << ia[i] << std::endl;
      }
      for( i = 0; i < NNZ; i++ )
      {
         mat_file << ja[i] << std::endl;
      }
      for( i = 0; i < NNZ; i++ )
      {
         mat_file << a_[i] << std::endl;
      }

      /* Right hand side. */
      if( rhs_vals )
         for( i = 0; i < N; i++ )
         {
            mat_file << rhs_vals[i] << std::endl;
         }
   }

   /* additional matrix format */
   if( getenv("IPOPT_WRITE_MAT_MTX") )
   {
      /* Write header */
      char mat_name[128];
      char mat_pref[32];

      Index i;
      Index j;

      if( getenv("IPOPT_WRITE_PREFIX") )
      {
         strcpy(mat_pref, getenv("IPOPT_WRITE_PREFIX"));
      }
      else
      {
         strcpy(mat_pref, "mat-ipopt");
      }

      Snprintf(mat_name, 127, "%s_%03d-%02d.mtx", mat_pref, iter_cnt, sol_cnt);

      // Open and write matrix file.
      std::ofstream mat_file(mat_name);
      mat_file << std::setprecision(std::numeric_limits<Number>::digits10 + 1);

      for( i = 0; i < N; i++ )
         for( j = ia[i]; j < ia[i + 1] - 1; j++ )
         {
            mat_file << ' ' << i + 1 << ' ' << ja[j - 1] << ' ' << a_[j - 1] << std::endl;
         }
   }
}

ESymSolverStatus PardisoSolverInterface::Factorization(
   const Index* ia,
   const Index* ja,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   DBG_START_METH("PardisoSolverInterface::Factorization", dbg_verbosity);

   // Call Pardiso to do the factorization
   Index PHASE;
   Index N = dim_;
   Index PERM;   // This should not be accessed by Pardiso
   Index NRHS = 0;
   Number B;  // This should not be accessed by Pardiso in factorization
   // phase
   Number X;  // This should not be accessed by Pardiso in factorization
   // phase
   Index ERROR;

   bool done = false;
   bool just_performed_symbolic_factorization = false;

   while( !done )
   {
      if( !have_symbolic_factorization_ )
      {
         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
         }
         PHASE = 11;

#ifdef PARDISO_MATCHING_PREPROCESS
         delete[] ia2;
         ia2 = NULL;

         delete[] ja2;
         ja2 = NULL;

         delete[] a2_;
         a2_ = NULL;

         delete[] perm2;
         perm2 = NULL;

         delete[] scale2;
         scale2 = NULL;

         ia2 = new Index[N + 1];
         ja2 = new Index[nonzeros_];
         a2_ = new Number[nonzeros_];
         perm2 = new Index[N];
         scale2 = new Number[N];
         Index* tmp2_ = new Index[N];

         smat_reordering_pardiso_wsmp(&N, ia, ja, a_, ia2, ja2, a2_, perm2, scale2, tmp2_, 0);

         delete[] tmp2_;

#endif

         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Calling Pardiso for symbolic factorization.\n");
         pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
#ifdef PARDISO_MATCHING_PREPROCESS
                 &PHASE, &N, a2_, ia2, ja2, &PERM,
#else
                 &PHASE, &N, a_, ia, ja, &PERM,
#endif
                 &NRHS, IPARM_, &MSGLVL_, &B, &X, &ERROR, DPARM_);
         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemSymbolicFactorization().End();
         }
         if( ERROR == -7 )
         {
            Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                           "Pardiso symbolic factorization returns ERROR = %" IPOPT_INDEX_FORMAT ".  Matrix is singular.\n", ERROR);
            return SYMSOLVER_SINGULAR;
         }
         else if( ERROR != 0 )
         {
            Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                           "Error in Pardiso during symbolic factorization phase.  ERROR = %" IPOPT_INDEX_FORMAT ".\n", ERROR);
            return SYMSOLVER_FATAL_ERROR;
         }
         have_symbolic_factorization_ = true;
         just_performed_symbolic_factorization = true;

         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Memory in KB required for the symbolic factorization  = %" IPOPT_INDEX_FORMAT ".\n", IPARM_[14]);
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Integer memory in KB required for the numerical factorization  = %" IPOPT_INDEX_FORMAT ".\n", IPARM_[15]);
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Double  memory in KB required for the numerical factorization  = %" IPOPT_INDEX_FORMAT ".\n", IPARM_[16]);
      }

      PHASE = 22;

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().Start();
      }
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "Calling Pardiso for factorization.\n");
      // Dump matrix to file, and count number of solution steps.
      if( HaveIpData() )
      {
         if( IpData().iter_count() != debug_last_iter_ )
         {
            debug_cnt_ = 0;
         }
         debug_last_iter_ = IpData().iter_count();
         debug_cnt_++;
      }
      else
      {
         debug_cnt_ = 0;
         debug_last_iter_ = 0;
      }

#ifdef PARDISO_MATCHING_PREPROCESS
      Index* tmp3_ = new Index[N];
      smat_reordering_pardiso_wsmp(&N, ia, ja, a_, ia2, ja2, a2_, perm2, scale2, tmp3_, 1);
      delete[] tmp3_;
#endif

      pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
#ifdef PARDISO_MATCHING_PREPROCESS
              &PHASE, &N, a2_, ia2, ja2, &PERM,
#else
              &PHASE, &N, a_, ia, ja, &PERM,
#endif
              &NRHS, IPARM_, &MSGLVL_, &B, &X, &ERROR, DPARM_);
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }

      if( ERROR == -7 )
      {
         Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                        "Pardiso factorization returns ERROR = %" IPOPT_INDEX_FORMAT ".  Matrix is singular.\n", ERROR);
         return SYMSOLVER_SINGULAR;
      }
      else if( ERROR == -4 )
      {
         // I think this means that the matrix is singular
         // OLAF said that this will never happen (ToDo)
         return SYMSOLVER_SINGULAR;
      }
      else if( ERROR != 0 )
      {
         Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                        "Error in Pardiso during factorization phase.  ERROR = %" IPOPT_INDEX_FORMAT ".\n", ERROR);
         return SYMSOLVER_FATAL_ERROR;
      }

      negevals_ = Max(IPARM_[22], numberOfNegEVals);
      if( IPARM_[13] != 0 )
      {
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "Number of perturbed pivots in factorization phase = %" IPOPT_INDEX_FORMAT ".\n", IPARM_[13]);
         if( !pardiso_redo_symbolic_fact_only_if_inertia_wrong_ || (negevals_ != numberOfNegEVals) )
         {
            if( HaveIpData() )
            {
               IpData().Append_info_string("Pn");
            }
            have_symbolic_factorization_ = false;
            // We assume now that if there was just a symbolic
            // factorization and we still have perturbed pivots, that
            // the system is actually singular, if
            // pardiso_repeated_perturbation_means_singular_ is true
            if( just_performed_symbolic_factorization )
            {
               if( pardiso_repeated_perturbation_means_singular_ )
               {
                  if( HaveIpData() )
                  {
                     IpData().Append_info_string("Ps");
                  }
                  return SYMSOLVER_SINGULAR;
               }
               else
               {
                  done = true;
               }
            }
            else
            {
               done = false;
            }
         }
         else
         {
            if( HaveIpData() )
            {
               IpData().Append_info_string("Pp");
            }
            done = true;
         }
      }
      else
      {
         done = true;
      }
   }

   DBG_ASSERT(IPARM_[21] + IPARM_[22] == dim_);

   // Check whether the number of negative eigenvalues matches the requested
   // count
   if( skip_inertia_check_ )
   {
      numberOfNegEVals = negevals_;
   }

   if( check_NegEVals && (numberOfNegEVals != negevals_) )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %" IPOPT_INDEX_FORMAT ", but we got %" IPOPT_INDEX_FORMAT ".\n", numberOfNegEVals, negevals_);
      return SYMSOLVER_WRONG_INERTIA;
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus PardisoSolverInterface::Solve(
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   Number*      rhs_vals
)
{
   DBG_START_METH("PardisoSolverInterface::Solve", dbg_verbosity);

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().Start();
   }
   // Call Pardiso to do the solve for the given right-hand sides
   Index PHASE = 33;
   Index N = dim_;
   Index PERM;   // This should not be accessed by Pardiso
   Index NRHS = nrhs;
   Number* X = new Number[nrhs * dim_];

   Number* ORIG_RHS = new Number[nrhs * dim_];
   Index ERROR;
   // Initialize solution with zero and save right hand side
   for( Index i = 0; i < N; i++ )
   {
      X[i] = 0.;
      ORIG_RHS[i] = rhs_vals[i];
   }

   // Dump matrix to file if requested
   Index iter_count = 0;
   if( HaveIpData() )
   {
      iter_count = IpData().iter_count();
   }

#ifdef PARDISO_MATCHING_PREPROCESS
   write_iajaa_matrix (N, ia2, ja2, a2_, rhs_vals, iter_count, debug_cnt_);
#else
   write_iajaa_matrix(N, ia, ja, a_, rhs_vals, iter_count, debug_cnt_);
#endif

   int attempts = 0;
   const int max_attempts = pardiso_iterative_ ? pardiso_max_droptol_corrections_ + 1 : 1;
   DBG_ASSERT(max_attempts > 0);

   while( attempts < max_attempts )
   {

#ifdef PARDISO_MATCHING_PREPROCESS
      for (Index i = 0; i < N; i++)
      {
         rhs_vals[perm2[i]] = scale2[i] * ORIG_RHS[ i ];
      }
      pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_,
              &PHASE, &N, a2_, ia2, ja2, &PERM,
              &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
              &ERROR, DPARM_);
      for (Index i = 0; i < N; i++)
      {
         X[i] = rhs_vals[ perm2[i]];
      }
      for (Index i = 0; i < N; i++)
      {
         rhs_vals[i] = scale2[i] * X[i];
      }

#else
      for( Index i = 0; i < N; i++ )
      {
         rhs_vals[i] = ORIG_RHS[i];
      }
      pardiso(PT_, &MAXFCT_, &MNUM_, &MTYPE_, &PHASE, &N, a_, ia, ja, &PERM, &NRHS, IPARM_, &MSGLVL_, rhs_vals, X,
              &ERROR, DPARM_);
#endif

      if( ERROR <= -100 && ERROR >= -102 )
      {
         Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                        "Iterative solver in Pardiso did not converge (ERROR = %" IPOPT_INDEX_FORMAT ")\n", ERROR);
         Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                        "  Decreasing drop tolerances from DPARM_[4] = %e and DPARM_[5] = %e\n", DPARM_[4], DPARM_[5]);
         PHASE = 23;
         DPARM_[4] /= 2.0;
         DPARM_[5] /= 2.0;
         Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                        "                               to DPARM_[4] = %e and DPARM_[5] = %e\n", DPARM_[4], DPARM_[5]);
         attempts++;
         ERROR = 0;
      }
      else
      {
         attempts = max_attempts;
         // TODO we could try again with some PARDISO parameters changed, i.e., enabling iterative refinement
      }
   }

   delete[] X;
   delete[] ORIG_RHS;

   if( IPARM_[6] != 0 )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of iterative refinement steps = %" IPOPT_INDEX_FORMAT ".\n", IPARM_[6]);
      if( HaveIpData() )
      {
         IpData().Append_info_string("Pi");
      }
   }

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemBackSolve().End();
   }
   if( ERROR != 0 )
   {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error in Pardiso during solve phase.  ERROR = %" IPOPT_INDEX_FORMAT ".\n", ERROR);
      return SYMSOLVER_FATAL_ERROR;
   }
   return SYMSOLVER_SUCCESS;
}

Index PardisoSolverInterface::NumberOfNegEVals() const
{
   DBG_START_METH("PardisoSolverInterface::NumberOfNegEVals", dbg_verbosity);
   DBG_ASSERT(negevals_ >= 0);
   return negevals_;
}

bool PardisoSolverInterface::IncreaseQuality()
{
   // At the moment, I don't see how we could tell Pardiso to do better
   // (maybe switch from IPARM[20]=1 to IPARM[20]=2?)
   return false;
}

} // namespace Ipopt
