
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <cuda_runtime.h>
#include "cudss.h"

#include "cuDSS_wrapper.h"

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
static cudssMatrix_t a_;
static cudssMatrix_t b_;
static cudssMatrix_t sol_;
// static const cudssMatrixFormat_t matFormat_ = CUDSS_MFORMAT_CSR;
static const cudssMatrixType_t matType_ = CUDSS_MTYPE_GENERAL;
static const cudssMatrixViewType_t matViewType_ = CUDSS_MVIEW_UPPER;
static const cudssIndexBase_t matIndex_ = CUDSS_BASE_ZERO;

void cuDSS_initialize() {
    // ADD ERROR CHECKING EVERYWHERE FOR CUDA AND CUDSS
    cudaStreamCreate(&stream_);
    status_ = cudssCreate(&handle_);
    status_ = cudssSetStream(handle_, stream_);
}

void cuDSS_terminate() {
    status_ = cudssMatrixDestroy(a_);
    status_ = cudssMatrixDestroy(b_);
    status_ = cudssMatrixDestroy(sol_);
    status_ = cudssDataDestroy(handle_, data_);
    status_ = cudssConfigDestroy(config_);
    status_ = cudssDestroy(handle_);
    cudaStreamSynchronize(stream_);
    cudaFree(aD_);
    cudaFree(ia_);
    cudaFree(ja_);
}

bool cuDSS_config_create_and_set(cuDSS_config_settings settings) {
    // Creating cuDSS solver configuration
    status_ = cudssConfigCreate(&config_);

    algReorder_ = static_cast<cudssAlgType_t>(settings.algReorder);
    algFactor_ = static_cast<cudssAlgType_t>(settings.algFactor);
    algPivotEps_ = static_cast<cudssAlgType_t>(settings.algPivotEps);
    algMatching_ = static_cast<cudssAlgType_t>(settings.algMatching);
    pivotType_ = static_cast<cudssPivotType_t>(settings.pivotType);

    useMatching_ = settings.useMatching;
    nIterSteps_ = settings.nIterSteps;
    pivotThr_ = settings.pivotThr;
    pivotEps_ = settings.pivotEps;
    maxLUnnz_ = settings.maxLUnnz;
    nThreads_ = settings.nThreads;
    ndNLevels_ = settings.ndNLevels;
    useSP_ = settings.useSP;
    schurMode_ = settings.schurMode;
    deterministic_ = settings.deterministic;

    status_ = cudssConfigSet(config_, CUDSS_CONFIG_REORDERING_ALG, &algReorder_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_FACTORIZATION_ALG, &algFactor_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_EPSILON_ALG, &algPivotEps_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_USE_MATCHING, &useMatching_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_MATCHING_ALG, &algMatching_, sizeof(cudssAlgType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_IR_N_STEPS, &nIterSteps_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_TYPE, &pivotType_, sizeof(cudssPivotType_t));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_THRESHOLD, &pivotThr_, sizeof(double));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_PIVOT_EPSILON, &pivotEps_, sizeof(double));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_MAX_LU_NNZ, &maxLUnnz_, sizeof(Index));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_HOST_NTHREADS, &nThreads_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_ND_NLEVELS, &ndNLevels_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_USE_SUPERPANELS, &useSP_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_SCHUR_MODE, &schurMode_, sizeof(int));
    status_ = cudssConfigSet(config_, CUDSS_CONFIG_DETERMINISTIC_MODE, &deterministic_, sizeof(int));

    // Creating cuDSS data container
    status_ = cudssDataCreate(handle_, &data_);

    if (status_ != CUDSS_STATUS_SUCCESS) return false;

    return true;

}

void cuDSS_initialize_structure(Index dim, Index nonzeros, const Index* ia, const Index* ja) {
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Storing the row and column indexes on Device
    ia_ = NULL;
    cudaMalloc(&ia_, (dim_ + 1) * sizeof(Index));
    cudaMemcpy(ia_, ia, (dim_ + 1) * sizeof(Index), cudaMemcpyHostToDevice);
    ja_ = NULL;
    cudaMalloc(&ja_, nonzeros_ * sizeof(Index));
    cudaMemcpy(ja_, ja, nonzeros_ * sizeof(Index), cudaMemcpyHostToDevice);

    // Storing the matrix elements on Device and Host
    Number* aH = NULL;
    aH = new Number[nonzeros_];
    aD_ = NULL;
    cudaMalloc(&aD_, nonzeros_ * sizeof(Number));
    cudaMemcpy(aD_, aH, nonzeros_ * sizeof(Number), cudaMemcpyHostToDevice);
    status_ = cudssMatrixCreateCsr(&a_, dim_, dim_, nonzeros_, ia_, NULL,
                                    ja_, aD_, CUDA_R_32I, CUDA_R_64F, matType_, 
                                    matViewType_, matIndex_);
    delete[] aH;
}

int cuDSS_reordering()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_REORDERING, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) return 4;
    return 0;
}

int cuDSS_symbolic_factorization()
{
    status_ = cudssExecute(handle_, CUDSS_PHASE_SYMBOLIC_FACTORIZATION, config_, data_, a_, sol_, b_);
    if (status_ != CUDSS_STATUS_SUCCESS) return 4;
    return 0;
}

Number* cuDSS_get_matrix_values()
{
    return aD_;
}
