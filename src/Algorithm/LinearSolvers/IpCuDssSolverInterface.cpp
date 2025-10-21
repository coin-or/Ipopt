// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpCuDssSolverInterface.hpp"

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

cuDSSSolverInterface::cuDSSSolverInterface()
{
    cudaStreamCreate(&stream_);
    status_ = cudssCreate(&handle_);
    status_ = cudssSetStream(handle_, stream_);
}

cuDSSSolverInterface::~cuDSSSolverInterface()
{
    status_ = cudssMatrixDestroy(A);
    status_ = cudssMatrixDestroy(b);
    status_ = cudssMatrixDestroy(x);
    status_ = cudssDataDestroy(handle_, data_);
    status_ = cudssConfigDestroy(config_);
    status_ = cudssDestroy(handle_);
    cudaStreamSynchronize(stream_);

    delete[] aH_;
    delete[] aD_;
}

void cuDSSSolverInterface::RegisterOptions(
    SmartPtr<RegisteredOptions> roptions
)
{
    roptions->AddStringOption4(
        "cuDSS_reordering_algorithm",
        "Algorithm for the reordering phase used by cuDSS",
        "CUDSS_ALG_DEFAULT",
        "CUDSS_ALG_DEFAULT", "Metis Based",
        "CUDSS_ALG_1", "Colamd based",
        "CUDSS_ALG_2", "Colamd based + special factorization",
        "CUDSS_ALG_3", "Amd based");
    roptions->AddStringOption2(
        "cuDSS_factorization_algorithm",
        "Algorithm for the factorization phase used by cuDSS",
        "CUDSS_ALG_DEFAULT",
        "CUDSS_ALG_DEFAULT", "Default",
        "CUDSS_ALG_1", "Modified default");
    roptions->AddStringOption2(
        "cuDSS_pivot_epsilon_algorithm",
        "Algorithm for the pivot epsilon calculation used by cuDSS",
        "CUDSS_ALG_DEFAULT",
        "CUDSS_ALG_DEFAULT", "Replace",
        "CUDSS_ALG_1", "Replace and Scale");
    roptions->AddBoolOption(
        "cuDSS_use_matching",
        "Flag for enabling/disabling matching.",
        false);
    roptions->AddStringOption6(
        "cuDSS_matching_algorithm",
        "Algorithm for matching used by cuDSS",
        "CUDSS_ALG_DEFAULT",
        "CUDSS_ALG_DEFAULT", "Default",
        "CUDSS_ALG_1", "",
        "CUDSS_ALG_2", "",
        "CUDSS_ALG_3", "",
        "CUDSS_ALG_4", "",
        "CUDSS_ALG_5", "");
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_number_iterative_steps",
        "Number of steps during the iterative refinement.",
        0,
        0);
    roptions->AddStringOption3(
        "cuDSS_pivot_type",
        "Type of pivoting used by cuDSS",
        "CUDSS_PIVOT_COL",
        "CUDSS_PIVOT_COL", "Default",
        "CUDSS_PIVOT_ROW", "",
        "CUDSS_PIVOT_NONE", "");
    roptions->AddLowerBoundedNumberOption(
        "cuDSS_pivoting_threshold",
        "Pivoting threshold.",
        0.0f,
        true,
        1.0f);
    roptions->AddLowerBoundedNumberOption(
        "cuDSS_pivoting_epsilon",
        "Pivoting epsilon.",
        0.0,
        true,
        #ifdef IPOPT_SINGLE
            1e-5
        #else
            1e-13
        #endif
        );
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_max_nnz_LU",
        "Upper limit on the number of nonzero entries in LU factors.",
        -1,
        -1);
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_nThreads",
        "Number of threads to be used by cuDSS in MT mode.",
        -1,
        -1);
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_min_NDLevels",
        "Minimum number of levels for the nested dissection reordering.",
        1,
        10);
    roptions->AddBoolOption(
        "cuDSS_use_superpanels",
        "Use superpanel optimization: 1 (default = enabled) or 0 (disabled).",
        true);
    roptions->AddBoolOption(
        "cuDSS_schur_mode",
        "Schur complement mode: 0 (default = disabled) or 1 (enabled).",
        false);
    roptions->AddBoolOption(
        "cuDSS_determ_mode",
        "Enable deterministic mode.",
        false);
}

bool cuDSSSolverInterface::InitializeImpl(
    const OptionsList &options,
    const std::string &prefix
)
{
    Index algReorder;
    options.GetEnumValue("cuDSS_reordering_algorithm", algReorder, prefix);
    algReorder_ = static_cast<cudssAlgType_t>(algReorder);
    Index algFactor;
    options.GetEnumValue("cuDSS_factorization_algorithm", algFactor, prefix);
    algFactor_ = static_cast<cudssAlgType_t>(algFactor);
    Index algPivotEps;
    options.GetEnumValue("cuDSS_pivot_epsilon_algorithm", algPivotEps, prefix);
    algPivotEps_ = static_cast<cudssAlgType_t>(algPivotEps);
    bool useMatching;
    options.GetBoolValue("cuDSS_use_matching", useMatching, prefix);
    useMatching_ = static_cast<int>(useMatching);
    Index algMatching;
    options.GetEnumValue("cuDSS_matching_algorithm", algMatching, prefix);
    algMatching_ = static_cast<cudssAlgType_t>(algMatching);
    Index nIterSteps;
    options.GetIntegerValue("cuDSS_number_iterative_steps", nIterSteps, prefix);
    nIterSteps_ = static_cast<int>(nIterSteps);
    Index pivotType;
    options.GetEnumValue("cuDSS_pivot_type", pivotType, prefix);
    pivotType_ = static_cast<cudssPivotType_t>(pivotType);
    Number pivotThr;
    options.GetNumericValue("cuDSS_pivoting_threshold", pivotThr, prefix);
    pivotThr_ = static_cast<double>(pivotThr);
    Number pivotEps;
    options.GetNumericValue("cuDSS_pivoting_epsilon", pivotEps, prefix);
    pivotEps_ = static_cast<double>(pivotEps);
    Index maxLUnnz;
    options.GetIntegerValue("cuDSS_max_nnz_LU", maxLUnnz, prefix);
    maxLUnnz_ = maxLUnnz;
    Index nThreads;
    options.GetIntegerValue("cuDSS_nThreads", nThreads, prefix);
    nThreads_ = static_cast<int>(nThreads);
    Index ndNLevels;
    options.GetIntegerValue("cuDSS_min_NDLevels", ndNLevels, prefix);
    ndNLevels_ = static_cast<int>(ndNLevels);
    bool useSP;
    options.GetBoolValue("cuDSS_use_superpanels", useSP, prefix);
    useSP_ = static_cast<int>(useSP);
    bool schurMode;
    options.GetBoolValue("cuDSS_schur_mode", schurMode, prefix);
    schurMode_ = static_cast<int>(schurMode);
    bool deterministic;
    options.GetBoolValue("cuDSS_determ_mode", deterministic, prefix);
    deterministic_ = static_cast<int>(deterministic);

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "cuDSS matrix ordering CUDSS_CONFIG_REORDERING_ALG: %d\n", algReorder_);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "cuDSS matrix ordering CUDSS_CONFIG_FACTORIZATION_ALG: %d\n", algFactor_);

    return true;
}

ESymSolverStatus cuDSSSolverInterface::InitializeStructure(
    Index dim, 
    Index nonzeros, 
    const Index *ia, 
    const Index *ja
)
{
    DBG_START_METH("cuDSSSolverInterface::InitializeStructure", dbg_verbosity);
    dim_ = dim;
    nonzeros_ = nonzeros;

    // Make space for storing the matrix elements on Device and Host
    delete[] aH_;
    delete[] aD_;
    aH_ = NULL;
    aD_ = NULL;
    aH_ = new Number[nonzeros_];
    aD_ = new Number[nonzeros_];

    // Do the symbolic factorization
    ESymSolverStatus retval = SymbolicFactorization(ia, ja);
    cudssExecute(handle_, CUDSS_PHASE_REORDERING, cudssConfig_t config, 
                cudssConfig_t data, cudssMatrix_t matrix, 
                cudssMatrix_t solution, cudssMatrix_t rhs);
    cudssExecute(handle_, CUDSS_PHASE_SYMBOLIC_FACTORIZATION, cudssConfig_t config, 
                cudssConfig_t data, cudssMatrix_t matrix, 
                cudssMatrix_t solution, cudssMatrix_t rhs);
    if( retval != SYMSOLVER_SUCCESS )
    {
        return retval;
    }

    return retval;
}

} // namespace Ipopt