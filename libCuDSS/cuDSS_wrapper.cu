// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cuda_runtime.h>
#include "cudss.h"

#include "cuDSS_wrapper.h"

// ADD SUPPORT FOR MGMN, MG, MT, HYBRID MODE and IPOPT_SINGLE
// FIX cudaMALLOC and FREE ERROR, add support for generic matrix (and transform from symm to general in Ipopt)

// Error checking as in: https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuDSS/simple/simple.cpp
#define CUDA_CALL_AND_CHECK(call, msg) \
    do { \
        cuda_error_ = call; \
        if (cuda_error_ != cudaSuccess) { \
            printf("Example FAILED: CUDA API returned error = %d, details: " #msg "\n", cuda_error_); \
        } \
    } while(0);


#define CUDSS_CALL_AND_CHECK(call, status, msg) \
    do { \
        status = call; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        } \
    } while(0);

static cudaError_t cuda_error_;
static cudssStatus_t status_;
static cudssHandle_t handle_;
static cudaStream_t stream_;
static cudssConfig_t config_;
static cudssData_t data_;
static cudssAlgType_t algReorder_ = CUDSS_ALG_DEFAULT;
static cudssAlgType_t algFactor_ = CUDSS_ALG_DEFAULT;
// static const cudssAlgType_t algSolve_ = CUDSS_ALG_DEFAULT;
static cudssAlgType_t algPivotEps_ = CUDSS_ALG_DEFAULT;
static int useMatching_ = 0;
static cudssAlgType_t algMatching_ = CUDSS_ALG_DEFAULT;
// static const int solveMode_ = 0;
static int nIterSteps_ = 0;
static cudssPivotType_t pivotType_ = CUDSS_PIVOT_COL;
static double pivotThr_ = 1.0f;
#ifdef CUDSS_SINGLE
static double pivotEps_ = 1e-5;
#else
static double pivotEps_ = 1e-13;
#endif
static Index maxLUnnz_ = -1;
// static const int memMode_ = 0;
// static const int useCUDAregMem_ = 1;
static int nThreads_ = -1;
// static const int hybridMode_ = 0;
static int ndNLevels_ = 10;
// static const int uBatchSize_ = 1;
// static const int uBatchIdx = -1;
static int useSP_ = 1;
// static const int nGPUs = 1;
// static const int* listGPUs_ = NULL;
static int schurMode_ = 0;
static int deterministic_ = 0;
static Index dim_;
static Index nonzeros_;
static Number* aD_;
static Index* ia_;
static Index* ja_;
static Number* bD_;
static Number* solD_;
static Number* aH_;
static Number* bH_;
static Number* solH_;
static cudssMatrix_t a_;
static cudssMatrix_t b_;
static cudssMatrix_t sol_;
// static const cudssMatrixFormat_t matFormat_ = CUDSS_MFORMAT_CSR;
static const cudssMatrixType_t matType_ = CUDSS_MTYPE_SYMMETRIC;
static const cudssMatrixViewType_t matViewType_ = CUDSS_MVIEW_UPPER;
static const cudssIndexBase_t matIndex_ = CUDSS_BASE_ZERO;

void cuDSS_initialize() {
    // ADD ERROR CHECKING EVERYWHERE FOR CUDA AND CUDSS
    CUDA_CALL_AND_CHECK(cudaStreamCreate(&stream_), "stream create");
    status_ = cudssCreate(&handle_);
    status_ = cudssSetStream(handle_, stream_);
}

void cuDSS_terminate() {
    // Checked with NVIDIA compute-sanitizer --tool memcheck ./test/hs071_c
    // It seems better not to deallocate anything on GPU in case of a double IpoptSolve execution
    // Otherwise segfault... 
    
    // status_ = cudssMatrixDestroy(a_);
    // status_ = cudssMatrixDestroy(b_);
    // status_ = cudssMatrixDestroy(sol_);
    // status_ = cudssDataDestroy(handle_, data_);
    // status_ = cudssConfigDestroy(config_);
    // status_ = cudssDestroy(handle_);
    // CUDA_CALL_AND_CHECK(cudaFree(aD_), "free aD_");
    // CUDA_CALL_AND_CHECK(cudaFree(bD_), "free bD_");
    // CUDA_CALL_AND_CHECK(cudaFree(solD_), "free solD_");
    // CUDA_CALL_AND_CHECK(cudaFree(ia_), "free ia_");
    // CUDA_CALL_AND_CHECK(cudaFree(ja_), "free ja_");
    delete[] aH_;
    delete[] bH_;
    delete[] solH_;
    // CUDA_CALL_AND_CHECK(cudaStreamSynchronize(stream_), "stream sync");
    // CUDA_CALL_AND_CHECK(cudaStreamDestroy(stream_), "stream destroy");
    // CUDA_CALL_AND_CHECK(cudaDeviceSynchronize(), "Device sync");
}

bool cuDSS_config_create_and_set(cuDSS_config_settings settings) {
    status_ = cudssConfigCreate(&config_);

    // Error checking as documented in https://docs.nvidia.com/cuda/cudss/types.html
    algReorder_ = static_cast<cudssAlgType_t>(settings.algReorder);
    if (algReorder_ >= CUDSS_ALG_4) {
        printf("CUDSS_ALG_4 and CUDSS_ALG_5 are invalid choices for CUDSS_CONFIG_REORDERING_ALG.\n");
        return false;
    }
    if ((algReorder_ == CUDSS_ALG_1) || (algReorder_ == CUDSS_ALG_2)) { // Matrix set as symmetric
        printf("CUDSS_ALG_1 and CUDSS_ALG_2 are only supported for general (non-symmetric or non-hermitian) matrices.\n");
        return false;
    }

    algFactor_ = static_cast<cudssAlgType_t>(settings.algFactor);
    if (algFactor_ >= CUDSS_ALG_2) {
        printf("CUDSS_ALG_2, CUDSS_ALG_3, CUDSS_ALG_4 and CUDSS_ALG_5 are invalid choices for CUDSS_CONFIG_FACTORIZATION_ALG.\n");
        return false;
    }

    algPivotEps_ = static_cast<cudssAlgType_t>(settings.algPivotEps);
    if (algPivotEps_ >= CUDSS_ALG_2) {
        printf("CUDSS_ALG_2, CUDSS_ALG_3, CUDSS_ALG_4 and CUDSS_ALG_5 are invalid choices for CUDSS_CONFIG_PIVOT_EPSILON_ALG.\n");
        return false;
    }
    if ((algPivotEps_ == CUDSS_ALG_1) && ((algReorder_ == CUDSS_ALG_1) || (algReorder_ == CUDSS_ALG_2))) {
        printf("CUDSS_ALG_1 for CUDSS_CONFIG_PIVOT_EPSILON_ALG is not supported when CUDSS_CONFIG_REORDERING_ALG is set to CUDSS_ALG_1 or CUDSS_ALG_2.\n");
        return false;
    }

    algMatching_ = static_cast<cudssAlgType_t>(settings.algMatching);
    if ((algMatching_ == CUDSS_ALG_1) && ((algReorder_ == CUDSS_ALG_1) || (algReorder_ == CUDSS_ALG_2))) {
        printf("Matching is not supported for CUDSS_ALG_1 and CUDSS_ALG_2 reordering algorithms (which use global pivoting to make the solution more accurate) or distributed matrices.\n");
        return false;
    }

    pivotType_ = static_cast<cudssPivotType_t>(settings.pivotType);

    useMatching_ = settings.useMatching;
    nIterSteps_ = settings.nIterSteps;
    // pivotThr_ = settings.pivotThr; // FOR NOW DISABLED, WILL BE REINTRODUCED WHEN ADDING USER OPTION
    // if (algReorder_ == CUDSS_ALG_DEFAULT) {
    //     printf("This parameter is only supported when reordering algorithm is set to CUDSS_ALG_1 or CUDSS_ALG_2.\n");
    //     return false;
    // }

    pivotEps_ = settings.pivotEps;
    maxLUnnz_ = settings.maxLUnnz;
    nThreads_ = settings.nThreads; // To check when implementing MT mode.
    ndNLevels_ = settings.ndNLevels;
    if (algReorder_ != CUDSS_ALG_DEFAULT) {
        printf("This setting only works when reordering algorithm is CUDSS_ALG_DEFAULT.\n");
        return false;
    }

    useSP_ = settings.useSP;
    schurMode_ = settings.schurMode; // To recheck when implementing MGMN, or MG mode.
    if ((algReorder_ != CUDSS_ALG_DEFAULT) || (algFactor_ == CUDSS_ALG_1) || useMatching_) {
        printf("Currently not supported when CUDSS_ALG_1 or CUDSS_ALG_2 is used for reordering, when MGMN mode or multi-GPU mode is used, or, when CUDSS_ALG_1 is used for the factorization. It is also not supported when a user permutation is set, for uniform and non-uniform batches, or when matching is enabled.\n");
        return false;
    }

    deterministic_ = settings.deterministic;
    if (deterministic_) printf("Currently the feature is supported only for single-gpu, single rhs and with hybrid memory mode (CUDSS_CONFIG_HYBRID_MODE) and hybrid execute mode (CUDSS_CONFIG_HYBRID_EXECUTE_MODE) disabled.\n");

    status_ = cudssConfigSet(config_, CUDSS_CONFIG_REORDERING_ALG, &algReorder_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_FACTORIZATION_ALG, &algFactor_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_EPSILON_ALG, &algPivotEps_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_USE_MATCHING, &useMatching_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_MATCHING_ALG, &algMatching_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_IR_N_STEPS, &nIterSteps_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_TYPE, &pivotType_, sizeof(cudssPivotType_t));
    //status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_THRESHOLD, &pivotThr_, sizeof(double));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_EPSILON, &pivotEps_, sizeof(double));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_MAX_LU_NNZ, &maxLUnnz_, sizeof(Index));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_HOST_NTHREADS, &nThreads_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_ND_NLEVELS, &ndNLevels_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_USE_SUPERPANELS, &useSP_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_SCHUR_MODE, &schurMode_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_DETERMINISTIC_MODE, &deterministic_, sizeof(int));

    status_ = cudssDataCreate(handle_, &data_);

    if (status_ != CUDSS_STATUS_SUCCESS) return false;

    return true;

}

void cuDSS_initialize_structure(Index dim, Index nonzeros, const Index* ia, const Index* ja) {
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Storing the row and column indexes on Device
    ia_ = NULL;
    CUDA_CALL_AND_CHECK(cudaMalloc(&ia_, (dim_ + 1) * sizeof(Index)), "cudaMalloc ia_");
    CUDA_CALL_AND_CHECK(cudaMemcpy(ia_, ia, (dim_ + 1) * sizeof(Index), cudaMemcpyHostToDevice), "cudaMemcpy ia_");
    ja_ = NULL;
    CUDA_CALL_AND_CHECK(cudaMalloc(&ja_, nonzeros_ * sizeof(Index)), "cudaMalloc ja_");
    CUDA_CALL_AND_CHECK(cudaMemcpy(ja_, ja, nonzeros_ * sizeof(Index), cudaMemcpyHostToDevice), "cudaMemcpy ja_");

    // Storing the matrix elements on Device and Host
    aH_ = NULL;
    aH_ = new Number[nonzeros_];
    aD_ = NULL;
    CUDA_CALL_AND_CHECK(cudaMalloc(&aD_, nonzeros_ * sizeof(Number)), "cudaMalloc aD_");
    CUDA_CALL_AND_CHECK(cudaMemcpy(aD_, aH_, nonzeros_ * sizeof(Number), cudaMemcpyHostToDevice), "cudaMalloc aD_");
    #ifdef CUDSS_SINGLE
    status_ = cudssMatrixCreateCsr( &a_, dim_, dim_, nonzeros_, ia_, NULL,
                                    ja_, aD_, CUDA_R_64I, CUDA_R_32F, matType_, 
                                    matViewType_, matIndex_);
    #else
    status_ = cudssMatrixCreateCsr( &a_, dim_, dim_, nonzeros_, ia_, NULL,
                                    ja_, aD_, CUDA_R_64I, CUDA_R_64F, matType_, 
                                    matViewType_, matIndex_);
    #endif

    // Storing the right hand side elements on Device and Host
    bH_ = NULL;
    bH_ = new Number[dim_];
    bD_ = NULL;
    CUDA_CALL_AND_CHECK(cudaMalloc(&bD_, dim_ * sizeof(Number)), "cudaMalloc bD_");
    CUDA_CALL_AND_CHECK(cudaMemcpy(bD_, bH_, dim_ * sizeof(Number), cudaMemcpyHostToDevice), "cudaMalloc bD_");
    #ifdef CUDSS_SINGLE
    status_ = cudssMatrixCreateDn(&b_, dim_, (Index)1, dim_, bD_, CUDA_R_32F, CUDSS_LAYOUT_COL_MAJOR);
    #else
    status_ = cudssMatrixCreateDn(&b_, dim_, (Index)1, dim_, bD_, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);
    #endif

    // Storing the solution elements on Device and Host
    solH_ = NULL;
    solH_ = new Number[dim_];
    solD_ = NULL;
    CUDA_CALL_AND_CHECK(cudaMalloc(&solD_, dim_ * sizeof(Number)), "cudaMalloc solD_");
    CUDA_CALL_AND_CHECK(cudaMemcpy(solD_, solH_, dim_ * sizeof(Number), cudaMemcpyHostToDevice), "cudaMalloc solD_");
    #ifdef CUDSS_SINGLE
    status_ = cudssMatrixCreateDn(&sol_, dim_, (Index)1, dim_, solD_, CUDA_R_32F, CUDSS_LAYOUT_COL_MAJOR);
    #else
    status_ = cudssMatrixCreateDn(&sol_, dim_, (Index)1, dim_, solD_, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR);
    #endif
}

int cuDSS_reordering()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_REORDERING, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) {
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: REORDERING\n", status_);
        return 4;
    }
    return 0;
}

int cuDSS_symbolic_factorization()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_SYMBOLIC_FACTORIZATION, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) {
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: SYM FACTORING\n", status_);
        return 4;
    }
    return 0;
}

Number* cuDSS_get_matrix_values()
{
    return aH_;
}

int cuDSS_factorization()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_FACTORIZATION, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) {
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: FACTORING\n", status_);
        return 4;
    }
    return 0;
}

int cuDSS_refactorization()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_REFACTORIZATION, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) {
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: REFACTOR\n", status_);
        return 4;
    }
    return 0;
}

int cuDSS_solve(Index nrhs, Number* rhs_vals)
{
    for (Index i = 0; i < nrhs; i++) {
        CUDA_CALL_AND_CHECK(cudaMemcpy(bD_, &rhs_vals[i * dim_], dim_ * sizeof(Number), cudaMemcpyHostToDevice), "cudaMemcpy bD_ to Device");
        status_ = cudssExecute(handle_, CUDSS_PHASE_SOLVE, config_, data_, a_, sol_, b_);
        if (status_ != CUDSS_STATUS_SUCCESS) {
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: SOLVE\n", status_);
            return 4;
        }
        CUDA_CALL_AND_CHECK(cudaMemcpy(&rhs_vals[i * dim_], solD_, dim_ * sizeof(Number), cudaMemcpyDeviceToHost), "cudaMemcpy sol_ to Host");
    }
    return 0;
}

Index cuDSS_get_inertia()
{   
    size_t sizeWritten;
    Index inertia[2];
    status_ = cudssDataGet(handle_, data_, CUDSS_DATA_INERTIA, &inertia, sizeof(inertia), &sizeWritten);
    return inertia[1];
}

bool cuDSS_update_matrix() 
{
    CUDA_CALL_AND_CHECK(cudaMemcpy(aD_, aH_, nonzeros_ * sizeof(Number), cudaMemcpyHostToDevice), "cudaMemcpy aH_ to Device");
    return true;
}