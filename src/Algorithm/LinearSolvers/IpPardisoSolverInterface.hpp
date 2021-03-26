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

#ifdef IPOPT_HAS_PARDISO_MKL
#define IPOPT_DECL_PARDISOINIT(x) void (x)( \
   void*         PT,    \
   const ipfint* MTYPE, \
   ipfint*       IPARM  \
)
#else
/* assuming PARDISO 4.0.0 or above */
#define IPOPT_DECL_PARDISOINIT(x) void (x)( \
   void*         PT,     \
   const ipfint* MTYPE,  \
   const ipfint* SOLVER, \
   ipfint*       IPARM,  \
   ipnumber*     DPARM,  \
   ipfint*       E       \
)
#endif

#define IPOPT_DECL_PARDISO(x) void (x)( \
   void**          PT,     \
   const ipfint*   MAXFCT, \
   const ipfint*   MNUM,   \
   const ipfint*   MTYPE,  \
   const ipfint*   PHASE,  \
   const ipfint*   N,      \
   const ipnumber* A,      \
   const ipfint*   IA,     \
   const ipfint*   JA,     \
   const ipfint*   PERM,   \
   const ipfint*   NRHS,   \
   ipfint*         IPARM,  \
   const ipfint*   MSGLVL, \
   ipnumber*       B,      \
   ipnumber*       X,      \
   ipfint*         E,      \
   ipnumber*       DPARM   \
)

#define IPOPT_DECL_SMAT_REORDERING_PARDISO_WSMP(x) void (x)( \
   const ipfint*   N,          \
   const ipfint*   ia,         \
   const ipfint*   ja,         \
   const ipnumber* a_,         \
   ipfint*         a2,         \
   ipfint*         ja2,        \
   ipnumber*       a2_,        \
   ipfint*         perm2,      \
   ipnumber*       scale2,     \
   ipfint*         tmp2_,      \
   ipfint          preprocess  \
)

namespace Ipopt
{

/** Interface to the linear solver Pardiso, derived from
 *  SparseSymLinearSolverInterface.
 */
class PardisoSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   PardisoSolverInterface(
      SmartPtr<LibraryLoader> pardisoloader_
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
   ipfint* ia2;
   ipfint* ja2;
   Number* a2_;
   ipfint* perm2;
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
   ipfint MAXFCT_;
   /** Actual matrix for the solution phase. Is always 1.*/
   ipfint MNUM_;
   /** Matrix type; real and symmetric indefinite.  Is always -2.*/
   ipfint MTYPE_;
   /** Parameter and info array for Pardiso. */
   ipfint* IPARM_;
   /** Parameter and info array for Pardiso. */
   Number* DPARM_;
   /** Message level. */
   ipfint MSGLVL_;
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
#endif
