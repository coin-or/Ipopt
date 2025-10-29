// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi  McLaren Automotive Ltd. 2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

/* some useful links:
* cuDSS documentation: https://docs.nvidia.com/cuda/cudss/
*/

#ifndef __IPCUDSSSOLVERINTERFACE_HPP__
#define __IPCUDSSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.hpp"

// External library definitions for NVIDIA cuDSS for Ipopt
#include "cuDSS_wrapper.h"

namespace Ipopt 
{

/** Interface to the linear solver cuDSS as distributed by NVIDIA, derived from SparseSymLinearSolverInterface.
 * The current implementation is developed for the single GPU case (both single and double precision is supported).
 * No MG (Multi-GPU) or MGMN (Multi-GPU-Multi-Node) support for now. Will be integrated in the future.
 * No Hybrid mode support for now. Will be integrated in the future.
 * @since 3.14.20
 */
class cuDSSSolverInterface : public SparseSymLinearSolverInterface 
{
  private:

      /** @name cuDSS structure containing configuration settings */
      ///@{
      /** cuDSS Configuration Settings */
      cuDSS_config_settings settings_;
      ///@}

      Index negevals_;

      bool initialized_;

      bool configured_;

      ESymSolverStatus Factorization(
        const Index* ia,
        const Index* ja,
        bool         check_NegEVals,
        Index        numberOfNegEVals
      );

      ESymSolverStatus Solve(
        const Index* ia,
        const Index* ja,
        Index        nrhs,
        Number*      rhs_vals
      );

  public:

      /** @name Constructor/Destructor */
      ///@{
      /** Constructor */
      cuDSSSolverInterface();

      /** Destructor */
      ~cuDSSSolverInterface();
      ///@}

      bool InitializeImpl(
        const OptionsList& options,
        const std::string& prefix
      );

      ESymSolverStatus InitializeStructure(
        Index        dim,
        Index        nonzeros,
        const Index* ia,
        const Index* ja
      );

      Number* GetValuesArrayPtr();

      ESymSolverStatus MultiSolve(
        bool         new_matrix,
        const Index* ia,
        const Index* ja,
        Index        nrhs,
        Number*      rhs_vals,
        bool         check_NegEVals,
        Index        numberOfNegEVals
      );

      Index NumberOfNegEVals() const;

      virtual bool IncreaseQuality()
      {
        return false;
      }

      virtual bool ProvidesInertia() const
      {
        return true;
      }

      EMatrixFormat MatrixFormat() const
      {
        return CSR_Format_0_Offset;
      }

      static void RegisterOptions(
        SmartPtr<RegisteredOptions> roptions
      );

};

}

#endif