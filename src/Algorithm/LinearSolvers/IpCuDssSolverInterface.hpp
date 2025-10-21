// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

/* some useful links:
* cuDSS documentation: https://docs.nvidia.com/cuda/cudss/
*/

#ifndef __IPCUDSSSOLVERINTERFACE_HPP__
#define __IPCUDSSSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

// NVIDIA CUDA and cuDSS libraries
#include "cudss.h"

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

      /** @name Information about cuDSS solver */
      ///@{
      /** cuDSS status */
      cudssStatus_t status_;

      /** cuDSS handle */
      cudssHandle_t handle_;

      /** CUDA/cuDSS stream */
      cudaStream_t stream_;

      /** cuDSS solver configuration */
      cudssConfig_t config_;

      /** cuDSS data */
      cudssData_t data_;
      ///@}

      /** @name Information about cuDSS algorithms */
      ///@{
      /** Algorithm for the Reordering phase */
      cudssAlgType_t algReorder_ = CUDSS_ALG_DEFAULT;

      /** Algorithm for the Factorization phase */
      cudssAlgType_t algFactor_ = CUDSS_ALG_DEFAULT;

      /** Algorithm for the Solve phase */
      const cudssAlgType_t algSolve_ = CUDSS_ALG_DEFAULT;

      /** Algorithm for the Pivot Epsilon calculation */
      cudssAlgType_t algPivotEps_ = CUDSS_ALG_DEFAULT;

      /** Use Matching Algorithm */
      int useMatching_ = 0;

      /** Algorithm for Matching calculation */
      cudssAlgType_t algMatching_ = CUDSS_ALG_DEFAULT;

      /** Potential modificator on the system matrix (e.g. transpose or conjugate transpose) */
      const int solveMode_ = 0;

      /** Number of steps during the iterative refinement */
      int nIterSteps_ = 0;

      /** Pivoting type definition */
      cudssPivotType_t pivotType_ = CUDSS_PIVOT_COL;

      /** Pivoting threshold */
      double pivotThr_ = 1.0f;

      /** Pivoting epsilon */
      #ifdef IPOPT_SINGLE
      double pivotEps_ = 1e-5;
      #else
      double pivotEps_ = 1e-13;
      #endif

      /** Upper limit on the number of nonzero entries in LU factors. */
      Index maxLUnnz_ = -1;

      /** Hybrid mode memory. */
      const int memMode_ = 0;

      /** Enable or disable usage of cudaHostRegister() by cuDSS hybrid memory mode. */
      const int useCUDAregMem_ = 1;

      /** Number of threads to be used by cuDSS in MT mode. */
      int nThreads_ = -1;

      /** Hybrid execute mode. */
      const int hybridMode_ = 0;

      /** Minimum number of levels for the nested dissection reordering. */
      int ndNLevels_ = 10;

      /** The number of matrices in a uniform batch of systems to be processed by cuDSS. */
      const int uBatchSize_ = 1;

      /** -1 or a 0-based index of matrix in a uniform batch which will be processed during factorization or solve phase. */
      const int uBatchIdx = -1;

      /** Use superpanel optimization */
      int useSP_ = 1;

      /** Number of devices (MG or MGMN) */
      const int nGPUs = 1;

      /** Device list (MG or MGMN) */
      const int* listGPUs_ = NULL;

      /** Schur complement mode. */
      int schurMode_ = 0;

      /** Deterministic mode. */
      int deterministic_ = 0;
      ///@}

      /** @name Information about the matrix */
      ///@{
      /** Number of rows and columns of the matrix */
      Index dim_;

      /** Number of nonzeros of the matrix in triplet representation. */
      Index nonzeros_;

      /** Array for storing the values of the matrix on the host. */
      Number* aH_;

      /** Array for storing the values of the matrix on the device. */
      Number* aD_;

      /** Matrix format definition */
      const cudssMatrixFormat_t matFormat_ = CUDSS_MFORMAT_CSR;

      /** Matrix type definition */
      const cudssMatrixType_t matType_ = CUDSS_MTYPE_GENERAL;

      /** Matrix view type definition */
      const cudssMatrixViewType_t matViewType_ = CUDSS_MVIEW_UPPER;

      /** Indexing base definition */
      const cudssIndexBase_t matType_ = CUDSS_BASE_ZERO;
      ///@}

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

      bool ProvidesDegeneracyDetection() const;

      ESymSolverStatus DetermineDependentRows(
        const Index*      /*ia*/,
        const Index*      /*ja*/,
        std::list<Index>& /*c_deps*/
      );

      static void RegisterOptions(
        SmartPtr<RegisteredOptions> roptions
      );

};

}

#endif