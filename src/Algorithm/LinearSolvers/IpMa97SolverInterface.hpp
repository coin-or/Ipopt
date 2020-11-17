// Copyright (C) 2012, The Science and Technology Facilities Council (STFC)
// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Jonathan Hogg                    STFC   2012-12-21
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPMA97SOLVERINTERFACE_HPP__
#define __IPMA97SOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
extern "C"
{
#include "hsl_ma97d.h"
}

namespace Ipopt
{

class Ma97SolverInterface: public SparseSymLinearSolverInterface
{
private:
   enum order_opts
   {
      ORDER_AUTO,
      ORDER_BEST,
      ORDER_AMD,
      ORDER_METIS,
      ORDER_MATCHED_AUTO,
      ORDER_MATCHED_AMD,
      ORDER_MATCHED_METIS
   };
   enum scale_opts
   {
      SWITCH_NEVER,
      SWITCH_AT_START,
      SWITCH_AT_START_REUSE,
      SWITCH_ON_DEMAND,
      SWITCH_ON_DEMAND_REUSE,
      SWITCH_NDELAY,
      SWITCH_NDELAY_REUSE,
      SWITCH_OD_ND,
      SWITCH_OD_ND_REUSE
   };

   int ndim_;            ///< Number of dimensions
   Number* val_;         ///< Storage for variables
   int numneg_;          ///< Number of negative pivots in last factorization
   int numdelay_;        ///< Number of delayed pivots last time we scaled
   void* akeep_;         ///< Stores pointer to factors (only understood Fortran code!)
   void* fkeep_;         ///< Stores pointer to factors (only understood Fortran code!)
   bool pivtol_changed_; ///< indicates if pivtol has been changed
   bool rescale_;        ///< Indicates if we should rescale next factorization
   Number* scaling_;     ///< Store scaling for reuse if doing dynamic scaling
   int fctidx_;          ///< Current factorization number to dump to

   /* Options */
   struct ma97_control control_;
   Number umax_;
   int ordering_;
   int scaling_type_;
   enum scale_opts switch_[3];
   int scaling_val_[3];
   int current_level_;
   bool dump_;

public:

   Ma97SolverInterface()
      : val_(NULL),
        numdelay_(0),
        akeep_(NULL),
        fkeep_(NULL),
        pivtol_changed_(false),
        rescale_(false),
        scaling_(NULL),
        fctidx_(0),
        scaling_type_(0),
        dump_(false)
   { }

   ~Ma97SolverInterface();

   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
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
   ///@}

   /** converts a scaling option name to its ma97 option number */
   static int ScaleNameToNum(
      const std::string& name
   );
};

} // namespace Ipopt

#endif
