// Copyright (C) 2011, Science and Technology Facilities Council
// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Jonathan Hogg                    STFC   2011-03-14
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPMA86SOLVERINTERFACE_HPP__
#define __IPMA86SOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpMa77SolverInterface.hpp"   // to get MC68 declaration macros
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

extern "C"
{
#ifdef IPOPT_SINGLE
#include "hsl_ma86s.h"
#else
#include "hsl_ma86d.h"
#endif
}

/// @since 3.14.0
#define IPOPT_DECL_MA86_DEFAULT_CONTROL(x) void (x)( \
   struct ma86_control* control \
)

/// @since 3.14.0
#define IPOPT_DECL_MA86_ANALYSE(x) void (x)( \
   const int                  n,       \
   const int                  ptr[],   \
   const int                  row[],   \
   int                        order[], \
   void**                     keep,    \
   const struct ma86_control* control, \
   struct ma86_info*          info     \
)

/// @since 3.14.0
#define IPOPT_DECL_MA86_FACTOR(x) void (x)( \
   const int                  n,       \
   const int                  ptr[],   \
   const int                  row[],   \
   const ipnumber             val[],   \
   const int                  order[], \
   void**                     keep,    \
   const struct ma86_control* control, \
   struct ma86_info*          info,    \
   const ipnumber             scale[]  \
)

/// @since 3.14.0
#define IPOPT_DECL_MA86_FACTOR_SOLVE(x) void (x)( \
   const int                  n,        \
   const int                  ptr[],    \
   const int                  row[],    \
   const ipnumber             val[],    \
   const int                  order[],  \
   void**                     keep,     \
   const struct ma86_control* control,  \
   struct ma86_info*          info,     \
   const int                  nrhs,     \
   const int                  ldx,      \
   ipnumber                   xx[],     \
   const ipnumber             scale[]   \
)

/// @since 3.14.0
#define IPOPT_DECL_MA86_SOLVE(x) void (x)( \
   const int                  job,    \
   const int                  nrhs,   \
   const int                  ldx,    \
   ipnumber*                  xx,     \
   const int                  order[],\
   void**                     keep,   \
   const struct ma86_control* control,\
   struct ma86_info*          info,   \
   const ipnumber             scale[] \
)

/// @since 3.14.0
#define IPOPT_DECL_MA86_FINALISE(x) void (x)( \
   void**                     keep,   \
   const struct ma86_control* control \
)

namespace Ipopt
{

class Ma86SolverInterface: public SparseSymLinearSolverInterface
{
private:
   enum order_opts
   {
      ORDER_AUTO,
      ORDER_AMD,
      ORDER_METIS
   };

   int ndim_;     ///< Number of dimensions
   Number* val_;  ///< Storage for variables
   int numneg_;   ///< Number of negative pivots in last factorization
   int* order_; ///< Fill reducing permutation
   void* keep_;   ///< Stores pointer to factors (only understood by Fortran code!)
   bool pivtol_changed_; ///< indicates if pivtol has been changed

   /* Options */
   struct ma86_control control_;
   Number umax_;
   int ordering_;

   /**@name MA68 and MC68 function pointers
    * @{
    */
   SmartPtr<LibraryLoader> hslloader;

   IPOPT_DECL_MA86_DEFAULT_CONTROL(*ma86_default_control);
   IPOPT_DECL_MA86_ANALYSE(*ma86_analyse);
   IPOPT_DECL_MA86_FACTOR(*ma86_factor);
   IPOPT_DECL_MA86_FACTOR_SOLVE(*ma86_factor_solve);
   IPOPT_DECL_MA86_SOLVE(*ma86_solve);
   IPOPT_DECL_MA86_FINALISE(*ma86_finalise);
   IPOPT_DECL_MC68_DEFAULT_CONTROL(*mc68_default_control);
   IPOPT_DECL_MC68_ORDER(*mc68_order);
   ///@}

public:

   Ma86SolverInterface(
      SmartPtr<LibraryLoader> hslloader_  ///< @since 3.14.0
   )  : val_(NULL),
      order_(NULL),
      keep_(NULL),
      pivtol_changed_(false),
      hslloader(hslloader_),
      ma86_default_control(NULL),
      ma86_analyse(NULL),
      ma86_factor(NULL),
      ma86_factor_solve(NULL),
      ma86_solve(NULL),
      ma86_finalise(NULL),
      mc68_default_control(NULL),
      mc68_order(NULL)
   { }

   ~Ma86SolverInterface();

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /// set MA86 and MC68 functions to use for every instantiation of this class
   /// @since 3.14.0
   static void SetFunctions(
      IPOPT_DECL_MA86_DEFAULT_CONTROL(*ma86_default_control),
      IPOPT_DECL_MA86_ANALYSE(*ma86_analyse),
      IPOPT_DECL_MA86_FACTOR(*ma86_factor),
      IPOPT_DECL_MA86_FACTOR_SOLVE(*ma86_factor_solve),
      IPOPT_DECL_MA86_SOLVE(*ma86_solve),
      IPOPT_DECL_MA86_FINALISE(*ma86_finalise),
      IPOPT_DECL_MC68_DEFAULT_CONTROL(*mc68_default_control),
      IPOPT_DECL_MC68_ORDER(*mc68_order)
   );

   bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** @name Methods for requesting solution of the linear system. */
   ///@{
   ESymSolverStatus InitializeStructure(
      Index        dim,
      Index        nonzeros,
      const Index* ia,
      const Index* ja
   );

   Number* GetValuesArrayPtr()
   {
      return val_;
   }

   ESymSolverStatus MultiSolve(
      bool         new_matrix,
      const Index* ia,
      const Index* ja,
      Index        nrhs,
      Number*      rhs_vals,
      bool         check_NegEVals,
      Index        numberOfNegEVals
   );

   Index NumberOfNegEVals() const
   {
      return numneg_;
   }
   ///@}

   //* @name Options of Linear solver */
   ///@{
   bool IncreaseQuality();

   bool ProvidesInertia() const
   {
      return true;
   }

   EMatrixFormat MatrixFormat() const
   {
      return CSR_Format_1_Offset;
   }
   ///@}

   /** @name Methods related to the detection of linearly dependent
    *  rows in a matrix */
   ///@{
   bool ProvidesDegeneracyDetection() const
   {
      return false;
   }

   ESymSolverStatus DetermineDependentRows(
      const Index*      /*ia*/,
      const Index*      /*ja*/,
      std::list<Index>& /*c_deps*/
   )
   {
      return SYMSOLVER_FATAL_ERROR;
   }
};

} // namespace Ipopt

#endif
