// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-17

#ifndef __IPMA27TSOLVERINTERFACE_HPP__
#define __IPMA27TSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

/// @since 3.14.0
#define IPOPT_DECL_MA27A(x) void (x)( \
   ipindex*       N,     \
   ipindex*       NZ,    \
   const ipindex* IRN,   \
   const ipindex* ICN,   \
   ipindex*       IW,    \
   ipindex*       LIW,   \
   ipindex*       IKEEP, \
   ipindex*       IW1,   \
   ipindex*       NSTEPS,\
   ipindex*       IFLAG, \
   ipindex*       ICNTL, \
   ipnumber*      CNTL,  \
   ipindex*       INFO,  \
   ipnumber*      OPS    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA27B(x) void (x)( \
   ipindex*       N,      \
   ipindex*       NZ,     \
   const ipindex* IRN,    \
   const ipindex* ICN,    \
   ipnumber*      A,      \
   ipindex*       LA,     \
   ipindex*       IW,     \
   ipindex*       LIW,    \
   ipindex*       IKEEP,  \
   ipindex*       NSTEPS, \
   ipindex*       MAXFRT, \
   ipindex*       IW1,    \
   ipindex*       ICNTL,  \
   ipnumber*      CNTL,   \
   ipindex*       INFO    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA27C(x) void (x)( \
   ipindex*   N,      \
   ipnumber*  A,      \
   ipindex*   LA,     \
   ipindex*   IW,     \
   ipindex*   LIW,    \
   ipnumber*  W,      \
   ipindex*   MAXFRT, \
   ipnumber*  RHS,    \
   ipindex*   IW1,    \
   ipindex*   NSTEPS, \
   ipindex*   ICNTL,  \
   ipnumber*  CNTL    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA27I(x) void (x)( \
   ipindex*  ICNTL, \
   ipnumber* CNTL   \
)

namespace Ipopt
{
/** Interface to the symmetric linear solver MA27, derived from
 *  SparseSymLinearSolverInterface.
 */
class Ma27TSolverInterface: public SparseSymLinearSolverInterface
{
public:
   /** @name Constructor/Destructor */
   ///@{
   /** Constructor */
   Ma27TSolverInterface(
      SmartPtr<LibraryLoader> hslloader_  ///< @since 3.14.0
   );

   /** Destructor */
   virtual ~Ma27TSolverInterface();
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
      const Index* airn,
      const Index* ajcn
   );

   virtual Number* GetValuesArrayPtr();

   virtual ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* airn,
      const Index* ajcn,
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
      return Triplet_Format;
   }
   ///@}

   ///@{
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );
   ///@}

   /// set MA27 functions to use for every instantiation of this class
   /// @since 3.14.0
   static void SetFunctions(
      IPOPT_DECL_MA27A(*ma27a),
      IPOPT_DECL_MA27B(*ma27b),
      IPOPT_DECL_MA27C(*ma27c),
      IPOPT_DECL_MA27I(*ma27i)
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
   Ma27TSolverInterface(
      const Ma27TSolverInterface&
   );

   /** Default Assignment Operator */
   void operator=(
      const Ma27TSolverInterface&
   );
   ///@}

   /**@name MA27 function pointers
    * @{
    */
   SmartPtr<LibraryLoader> hslloader;

   IPOPT_DECL_MA27A(*ma27a);
   IPOPT_DECL_MA27B(*ma27b);
   IPOPT_DECL_MA27C(*ma27c);
   IPOPT_DECL_MA27I(*ma27i);
   ///@}

   /** @name Information about the matrix */
   ///@{
   /** Number of rows and columns of the matrix */
   Index dim_;

   /** Number of nonzeros of the matrix */
   Index nonzeros_;
   ///@}

   /** @name Information about most recent factorization/solve */
   ///@{
   /** Number of negative eigenvalues */
   Index negevals_;
   ///@}

   /** @name Initialization flags */
   ///@{
   /** Flag indicating if internal data is initialized.
    *
    *  For initialization, this object needs to have seen a matrix.
    */
   bool initialized_;
   /** Flag indicating if the matrix has to be refactorized because
    *  the pivot tolerance has been changed.
    */
   bool pivtol_changed_;
   /** Flag that is true if we just requested the values of the
    *  matrix again (SYMSOLVER_CALL_AGAIN) and have to factorize
    *  again.
    */
   bool refactorize_;
   ///@}

   /** @name Solver specific data/options */
   ///@{
   /** Pivot tolerance */
   Number pivtol_;

   /** Maximal pivot tolerance */
   Number pivtolmax_;

   /** Factor for estimating initial value of liw */
   Number liw_init_factor_;
   /** Factor for estimating initial value of la */
   Number la_init_factor_;
   /** Factor for increaseing memory */
   Number meminc_factor_;
   /** Flag indicating whether the TNLP with identical structure has
    *  already been solved before.
    */
   bool warm_start_same_structure_;
   /** Flag indicating if the inertia is always assumed to be correct. */
   bool skip_inertia_check_;
   /** Flag indicating if MA27 should continue if a singular matrix
    * is detected, but right hands sides are still accepted.
    */
   bool ignore_singularity_;
   ///@}

   /** @name Data for the linear solver.
    * Storing factorization and other solver specific data structure.
    */
   ///@{
   /** integer control values */
   Index icntl_[30];
   /** real control values */
   Number cntl_[5];

   /** length of integer work space */
   Index liw_;
   /** integer work space */
   Index* iw_;

   /** MA27's IKEEP */
   Index* ikeep_;
   /** MA27's NSTEPS */
   Index nsteps_;
   /** MA27's MAXFRT */
   Index maxfrt_;

   /** length LA of A */
   Index la_;
   /** factor A of matrix */
   Number* a_;

   /** flag indicating that la should be increased before next factorization */
   bool la_increase_;
   /** flag indicating that liw should be increased before next factorization */
   bool liw_increase_;
   ///@}

   /** @name Internal functions */
   ///@{
   /** Call MA27AX and reserve memory for MA27 data.
    *
    *  Reserve memory for iw_ and ikeep_, call MA27AD to perform
    *  symbolic manipulations, and reserve all the remaining data memory
    */
   ESymSolverStatus SymbolicFactorization(
      const Index* airn,
      const Index* ajcn
   );

   /** Call MA27BX to factorize the Matrix.
    *
    *  It is assumed that the first nonzeros_ element of a_ contain the values
    *  of the matrix to be factorized.
    */
   ESymSolverStatus Factorization(
      const Index* airn,
      const Index* ajcn,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   /** Call MA27CX to do the backsolve. */
   ESymSolverStatus Backsolve(
      Index   nrhs,
      Number* rhs_vals
   );
   ///@}
};

} // namespace Ipopt
#endif
