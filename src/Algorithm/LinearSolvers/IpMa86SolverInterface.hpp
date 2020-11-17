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
extern "C"
{
#include "hsl_ma86d.h"
}

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
   Index* order_; ///< Fill reducing permutation
   void* keep_;   ///< Stores pointer to factors (only understood by Fortran code!)
   bool pivtol_changed_; ///< indicates if pivtol has been changed

   /* Options */
   struct ma86_control control_;
   Number umax_;
   int ordering_;

public:

   Ma86SolverInterface()
      : val_(NULL),
        order_(NULL),
        keep_(NULL),
        pivtol_changed_(false)
   { }

   ~Ma86SolverInterface();

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
};

} // namespace Ipopt

#endif
