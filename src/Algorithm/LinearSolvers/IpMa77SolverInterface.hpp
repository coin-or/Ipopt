// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Jonathan Hogg                    STFC   2013-30-05
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPMA77SOLVERINTERFACE_HPP__
#define __IPMA77SOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

extern "C"
{
#ifdef IPOPT_SINGLE
#include "hsl_ma77s.h"
#else
#include "hsl_ma77d.h"
#endif
#include "hsl_mc68i.h"
}

/// @since 3.14.0
#define IPOPT_DECL_MA77_DEFAULT_CONTROL(x) void (x)( \
   struct ma77_control* control \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_OPEN_NELT(x) void (x)( \
   const int                    n,      \
   const char*                  fname1, \
   const char*                  fname2, \
   const char*                  fname3, \
   const char*                  fname4, \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info,   \
   const int                    nelt    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_OPEN(x) void (x)( \
   const int                    n,      \
   const char*                  fname1, \
   const char*                  fname2, \
   const char*                  fname3, \
   const char*                  fname4, \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_INPUT_VARS(x) void (x)( \
   const int                    idx,    \
   const int                    nvar,   \
   const int                    list[], \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_INPUT_REALS(x) void (x)( \
   const int                    idx,    \
   const int                    length, \
   const ipnumber               reals[],\
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_ANALYSE(x) void (x)( \
   const int                    order[], \
   void**                       keep,    \
   const struct ma77_control*   control, \
   struct ma77_info*            info     \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_FACTOR(x) void (x)( \
   const int                    posdef, \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info,   \
   const ipnumber*              scale   \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_FACTOR_SOLVE(x) void (x)( \
   const int                    posdef, \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info,   \
   const ipnumber*              scale,  \
   const int                    nrhs,   \
   const int                    lx,     \
   ipnumber                     rhs[]   \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_SOLVE(x) void (x)( \
   const int                    job,    \
   const int                    nrhs,   \
   const int                    lx,     \
   ipnumber                     xx[],   \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info,   \
   const ipnumber*              scale   \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_RESID(x) void (x)( \
   const int                    nrhs,    \
   const int                    lx,      \
   const ipnumber               xx[],    \
   const int                    lresid,  \
   ipnumber                     resid[], \
   void**                       keep,    \
   const struct ma77_control*   control, \
   struct ma77_info*            info,    \
   ipnumber*                    anorm_bnd\
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_SCALE(x) void (x)( \
   ipnumber                     scale[], \
   void**                       keep,    \
   const struct ma77_control*   control, \
   struct ma77_info*            info,    \
   ipnumber*                    anorm    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_ENQUIRE_POSDEF(x) void (x)( \
   ipnumber                     d[],    \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info    \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_ENQUIRE_INDEF(x) void (x)( \
   int                          piv_order[], \
   ipnumber                     d[],         \
   void**                       keep,        \
   const struct ma77_control*   control,     \
   struct ma77_info*            info         \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_ALTER(x) void (x)( \
   const ipnumber               d[],     \
   void**                       keep,    \
   const struct ma77_control*   control, \
   struct ma77_info*            info     \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_RESTART(x) void (x)( \
   const char*                  restart_file, \
   const char*                  fname1,       \
   const char*                  fname2,       \
   const char*                  fname3,       \
   const char*                  fname4,       \
   void**                       keep,         \
   const struct ma77_control*   control,      \
   struct ma77_info*            info          \
)

/// @since 3.14.0
#define IPOPT_DECL_MA77_FINALISE(x) void (x)( \
   void**                       keep,   \
   const struct ma77_control*   control,\
   struct ma77_info*            info    \
)

/// @since 3.14.0
#define IPOPT_DECL_MC68_DEFAULT_CONTROL(x) void (x)( \
   struct mc68_control* control \
)

/// @since 3.14.0
#define IPOPT_DECL_MC68_ORDER(x) void (x)( \
   int                          ord,    \
   int                          n,      \
   const int                    ptr[],  \
   const int                    row[],  \
   int                          perm[], \
   const struct mc68_control*   control,\
   struct mc68_info*            info    \
)

namespace Ipopt
{

class Ma77SolverInterface: public SparseSymLinearSolverInterface
{
private:
   enum order_opts
   {
      ORDER_AMD,
      ORDER_METIS
   };

   int ndim_;    ///< Number of dimensions
   Number* val_; ///< Storage for variables
   int numneg_;  ///< Number of negative pivots in last factorization
   void* keep_;  ///< Stores pointer to factors (only understood by Fortran code!)
   bool pivtol_changed_; ///< indicates if pivtol has been changed

   /* Options */
   struct ma77_control control_;
   Number umax_;
   int ordering_;

   /**@name MA77 and MC68 function pointers
    * @{
    */
   SmartPtr<LibraryLoader> hslloader;

   IPOPT_DECL_MA77_DEFAULT_CONTROL(*ma77_default_control);
   IPOPT_DECL_MA77_OPEN_NELT(*ma77_open_nelt);
   IPOPT_DECL_MA77_OPEN(*ma77_open);
   IPOPT_DECL_MA77_INPUT_VARS(*ma77_input_vars);
   IPOPT_DECL_MA77_INPUT_REALS(*ma77_input_reals);
   IPOPT_DECL_MA77_ANALYSE(*ma77_analyse);
   IPOPT_DECL_MA77_FACTOR(*ma77_factor);
   IPOPT_DECL_MA77_FACTOR_SOLVE(*ma77_factor_solve);
   IPOPT_DECL_MA77_SOLVE(*ma77_solve);
   IPOPT_DECL_MA77_RESID(*ma77_resid);
   IPOPT_DECL_MA77_SCALE(*ma77_scale);
   IPOPT_DECL_MA77_ENQUIRE_POSDEF(*ma77_enquire_posdef);
   IPOPT_DECL_MA77_ENQUIRE_INDEF(*ma77_enquire_indef);
   IPOPT_DECL_MA77_ALTER(*ma77_alter);
   IPOPT_DECL_MA77_RESTART(*ma77_restart);
   IPOPT_DECL_MA77_FINALISE(*ma77_finalise);
   IPOPT_DECL_MC68_DEFAULT_CONTROL(*mc68_default_control);
   IPOPT_DECL_MC68_ORDER(*mc68_order);
   ///@}

public:

   Ma77SolverInterface(
      SmartPtr<LibraryLoader> hslloader_  ///< @since 3.14.0
   )  : val_(NULL),
      keep_(NULL),
      pivtol_changed_(false),
      hslloader(hslloader_),
      ma77_default_control(NULL),
      ma77_open_nelt(NULL),
      ma77_open(NULL),
      ma77_input_vars(NULL),
      ma77_input_reals(NULL),
      ma77_analyse(NULL),
      ma77_factor(NULL),
      ma77_factor_solve(NULL),
      ma77_solve(NULL),
      ma77_resid(NULL),
      ma77_scale(NULL),
      ma77_enquire_posdef(NULL),
      ma77_enquire_indef(NULL),
      ma77_alter(NULL),
      ma77_restart(NULL),
      ma77_finalise(NULL),
      mc68_default_control(NULL),
      mc68_order(NULL)
   { }

   ~Ma77SolverInterface();

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

   /// set MA77 and MC68 functions to use for every instantiation of this class
   /// @since 3.14.0
   static void SetFunctions(
      IPOPT_DECL_MA77_DEFAULT_CONTROL(*ma77_default_control),
      IPOPT_DECL_MA77_OPEN_NELT(*ma77_open_nelt),
      IPOPT_DECL_MA77_OPEN(*ma77_open),
      IPOPT_DECL_MA77_INPUT_VARS(*ma77_input_vars),
      IPOPT_DECL_MA77_INPUT_REALS(*ma77_input_reals),
      IPOPT_DECL_MA77_ANALYSE(*ma77_analyse),
      IPOPT_DECL_MA77_FACTOR(*ma77_factor),
      IPOPT_DECL_MA77_FACTOR_SOLVE(*ma77_factor_solve),
      IPOPT_DECL_MA77_SOLVE(*ma77_solve),
      IPOPT_DECL_MA77_RESID(*ma77_resid),
      IPOPT_DECL_MA77_SCALE(*ma77_scale),
      IPOPT_DECL_MA77_ENQUIRE_POSDEF(*ma77_enquire_posdef),
      IPOPT_DECL_MA77_ENQUIRE_INDEF(*ma77_enquire_indef),
      IPOPT_DECL_MA77_ALTER(*ma77_alter),
      IPOPT_DECL_MA77_RESTART(*ma77_restart),
      IPOPT_DECL_MA77_FINALISE(*ma77_finalise),
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
      return CSR_Full_Format_1_Offset;
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
