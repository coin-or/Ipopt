// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#ifndef __IPPARDISOSOLVERINTERFACE_HPP__
#define __IPPARDISOSOLVERINTERFACE_HPP__

#include "IpoptConfig.h"
#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

//#define PARDISO_MATCHING_PREPROCESS

/* assuming PARDISO 4.0.0 or above */
/// @since 3.14.0
#define IPOPT_DECL_PARDISOINIT(x) void (x)( \
   void*          PT,     \
   const ipindex* MTYPE,  \
   const ipindex* SOLVER, \
   ipindex*       IPARM,  \
   ipnumber*      DPARM,  \
   ipindex*       E       \
)

/// @since 3.14.0
#define IPOPT_DECL_PARDISO(x) void (x)( \
   void**           PT,     \
   const ipindex*   MAXFCT, \
   const ipindex*   MNUM,   \
   const ipindex*   MTYPE,  \
   const ipindex*   PHASE,  \
   const ipindex*   N,      \
   const ipnumber*  A,      \
   const ipindex*   IA,     \
   const ipindex*   JA,     \
   const ipindex*   PERM,   \
   const ipindex*   NRHS,   \
   ipindex*         IPARM,  \
   const ipindex*   MSGLVL, \
   ipnumber*        B,      \
   ipnumber*        X,      \
   ipindex*         E,      \
   ipnumber*        DPARM   \
)

/// @since 3.14.0
#define IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(x) void (x)( \
   const ipindex*   N,          \
   const ipindex*   ia,         \
   const ipindex*   ja,         \
   const ipnumber*  a_,         \
   ipindex*         a2,         \
   ipindex*         ja2,        \
   ipnumber*        a2_,        \
   ipindex*         perm2,      \
   ipnumber*        scale2,     \
   ipindex*         tmp2_,      \
   ipindex          preprocess  \
)

namespace Ipopt
{

/** Interface to the linear solver Pardiso as distributed by pardiso-project.org, derived from
 *  SparseSymLinearSolverInterface.
 */
class PardisoSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   PardisoSolverInterface(
      SmartPtr<LibraryLoader> pardisoloader_  ///< @since 3.14.0
   );

   /** Destructor */
   virtual ~PardisoSolverInterface();
   ///@}

   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** @name Methods for requesting solution of the linear system. */
   ///@{
   virtual ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* ia,
      const Index* ja
   );

   virtual Number* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      Number*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals);

   virtual Index NumberOfNegEVals() const;
   ///@}

   //* @name Options of Linear solver */
   ///@{
   virtual bool IncreaseQuality();

   virtual bool ProvidesInertia() const
   {
      return true;
   }

   EMatrixFormat MatrixFormat() const
   {
      return CSR_Format_1_Offset;
   }
   ///@}

   ///@{
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );
   ///@}

   /// set Pardiso functions to use for every instantiation of this class
   ///
   /// unless PARDISO_MATCHING_PREPROCESS has been defined, @arg smat_reordering_pardiso_wsmp is ignored
   /// @since 3.14.0
   static void SetFunctions(
      IPOPT_DECL_PARDISOINIT(*pardisoinit),
      IPOPT_DECL_PARDISO(*pardiso),
      bool isparallel,
      IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*smat_reordering_pardiso_wsmp)
#ifndef PARDISO_MATCHING_PREPROCESS
      = NULL
#endif
   );

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Copy Constructor */
   PardisoSolverInterface(
      const PardisoSolverInterface&);

   /** Default Assignment Operator */
   void operator=(
      const PardisoSolverInterface&);
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Number of nonzeros of the matrix in triplet representation. */
   Index nonzeros_;

   /** Array for storing the values of the matrix. */
   Number* a_;
   ///@}

#ifdef PARDISO_MATCHING_PREPROCESS
   /** Array for storing the values of a second matrix that has been already reordered. */
   Index* ia2;
   Index* ja2;
   Number* a2_;
   Index* perm2;
   Number* scale2;
#endif

   /** @name Information about most recent factorization/solve */
   ///@{
   /** Number of negative eigenvalues */
   Index negevals_;
   ///@}

   /** @name Solver specific options */
   ///@{
   /** Type for matching strategies */
   enum PardisoMatchingStrategy
   {
      COMPLETE,
      COMPLETE2x2,
      CONSTRAINT
   };
   /** Option that controls the matching strategy. */
   PardisoMatchingStrategy match_strat_;
   /** Flag indicating if symbolic factorization has already been performed. */
   bool have_symbolic_factorization_;
   /** Flag indicating whether the symbolic factorization should only
    *  be done after perturbed elements, if the inertia was wrong
    */
   bool pardiso_redo_symbolic_fact_only_if_inertia_wrong_;
   /** Flag indicating whether repeated perturbed elements even after
    *  a new symbolic factorization should be interpreted as a
    *  singular matrix
    */
   bool pardiso_repeated_perturbation_means_singular_;
   /** Flag indicating if the inertia is always assumed to be
    *  correct.
    */
   bool skip_inertia_check_;
   /** Flag indicating whether we are using the iterative solver in Pardiso. */
   bool pardiso_iterative_;
   /** Maximal number of decreases of drop tolerance during one solve. */
   Index pardiso_max_droptol_corrections_;
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *  For initialization, this object needs to have seen a matrix.
    */
   bool initialized_;
   ///@}

   /** @name Solver specific information */
   ///@{
   /** Internal data address pointers. */
   void** PT_;
   /** Maximal number of factors with identical nonzero
    *  structure. Here, we only store one factorization. Is always 1.
    */
   Index MAXFCT_;
   /** Actual matrix for the solution phase. Is always 1.*/
   Index MNUM_;
   /** Matrix type; real and symmetric indefinite.  Is always -2.*/
   Index MTYPE_;
   /** Parameter and info array for Pardiso. */
   Index* IPARM_;
   /** Parameter and info array for Pardiso. */
   Number* DPARM_;
   /** Message level. */
   Index MSGLVL_;
   ///@}

   /**@name Some counters for debugging */
   ///@{
   Index debug_last_iter_;
   Index debug_cnt_;
   ///@}

   /**@name PARDISO function pointers
    * @{
    */
   SmartPtr<LibraryLoader> pardisoloader;

   IPOPT_DECL_PARDISOINIT(*pardisoinit);
   IPOPT_DECL_PARDISO(*pardiso);
#ifdef PARDISO_MATCHING_PREPROCESS
   IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(*smat_reordering_pardiso_wsmp);
#endif

   bool pardiso_exist_parallel;
   /**@} */

   /** @name Internal functions */
   ///@{
   /** Call Pardiso to do the analysis phase. */
   ESymSolverStatus SymbolicFactorization(
      const Index* ia,
      const Index* ja
   );

   /** Call Pardiso to factorize the Matrix. */
   ESymSolverStatus Factorization(
      const Index* ia,
      const Index* ja,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   /** Call Pardiso to do the Solve. */
   ESymSolverStatus Solve(
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      Number*      rhs_vals
   );
   ///@}
};

} // namespace Ipopt

#endif  /* __IPPARDISOSOLVERINTERFACE_HPP__ */
