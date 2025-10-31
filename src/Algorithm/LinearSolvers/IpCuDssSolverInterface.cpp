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

#define IPOPT_CUDSS_DEBUG 0

cuDSSSolverInterface::cuDSSSolverInterface() : 
negevals_(-1),
initialized_(false),
configured_(false)
{
    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_initialize START\n");
    #endif

    cuDSS_initialize();

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_initialize END\n");
    #endif
}

cuDSSSolverInterface::~cuDSSSolverInterface()
{
    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_terminate START\n");
    #endif

    cuDSS_terminate();

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_terminate END\n");
    #endif
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
    const std::string &prefix)
{
    Index algReorder;
    options.GetEnumValue("cuDSS_reordering_algorithm", algReorder, prefix);
    settings_.algReorder = static_cast<int>(algReorder);
    Index algFactor;
    options.GetEnumValue("cuDSS_factorization_algorithm", algFactor, prefix);
    settings_.algFactor = static_cast<int>(algFactor);
    Index algPivotEps;
    options.GetEnumValue("cuDSS_pivot_epsilon_algorithm", algPivotEps, prefix);
    settings_.algPivotEps = static_cast<int>(algPivotEps);
    bool useMatching;
    options.GetBoolValue("cuDSS_use_matching", useMatching, prefix);
    settings_.useMatching = static_cast<int>(useMatching);
    Index algMatching;
    options.GetEnumValue("cuDSS_matching_algorithm", algMatching, prefix);
    settings_.algMatching = static_cast<int>(algMatching);
    Index nIterSteps;
    options.GetIntegerValue("cuDSS_number_iterative_steps", nIterSteps, prefix);
    settings_.nIterSteps = static_cast<int>(nIterSteps);
    Index pivotType;
    options.GetEnumValue("cuDSS_pivot_type", pivotType, prefix);
    settings_.pivotType = static_cast<int>(pivotType);
    Number pivotThr;
    options.GetNumericValue("cuDSS_pivoting_threshold", pivotThr, prefix);
    settings_.pivotThr = static_cast<double>(pivotThr);
    Number pivotEps;
    options.GetNumericValue("cuDSS_pivoting_epsilon", pivotEps, prefix);
    settings_.pivotEps = static_cast<double>(pivotEps);
    Index maxLUnnz;
    options.GetIntegerValue("cuDSS_max_nnz_LU", maxLUnnz, prefix);
    settings_.maxLUnnz = maxLUnnz;
    Index nThreads;
    options.GetIntegerValue("cuDSS_nThreads", nThreads, prefix);
    settings_.nThreads = static_cast<int>(nThreads);
    Index ndNLevels;
    options.GetIntegerValue("cuDSS_min_NDLevels", ndNLevels, prefix);
    settings_.ndNLevels = static_cast<int>(ndNLevels);
    bool useSP;
    options.GetBoolValue("cuDSS_use_superpanels", useSP, prefix);
    settings_.useSP = static_cast<int>(useSP);
    bool schurMode;
    options.GetBoolValue("cuDSS_schur_mode", schurMode, prefix);
    settings_.schurMode = static_cast<int>(schurMode);
    bool deterministic;
    options.GetBoolValue("cuDSS_determ_mode", deterministic, prefix);
    settings_.deterministic = static_cast<int>(deterministic);

    if (configured_) return true;

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_config_create_and_set START\n");
    #endif

    bool status = cuDSS_config_create_and_set(settings_);
    configured_ = true;

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_config_create_and_set END\n");
    #endif

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "cuDSS matrix ordering CUDSS_CONFIG_REORDERING_ALG: %d\n", settings_.algReorder);
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "cuDSS matrix ordering CUDSS_CONFIG_FACTORIZATION_ALG: %d\n", settings_.algFactor);

    return status;
}

ESymSolverStatus cuDSSSolverInterface::InitializeStructure(
    Index dim, 
    Index nonzeros, 
    const Index *ia, 
    const Index *ja
)
{
    DBG_START_METH("cuDSSSolverInterface::InitializeStructure", dbg_verbosity);
    
    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_initialize_structure START\n");
    #endif

    if (!initialized_) {
        cuDSS_initialize_structure(dim, nonzeros, ia, ja);

        #if IPOPT_CUDSS_DEBUG == 1
        printf("cuDSS_initialize_structure END\n");
        
        printf("cuDSS_symbolic_factorization START\n");
        #endif

        if (HaveIpData()) IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "cuDSS Symbolic Factorization.\n");
        if ( static_cast<ESymSolverStatus>(cuDSS_reordering()) != SYMSOLVER_SUCCESS ) return SYMSOLVER_FATAL_ERROR;
        if ( static_cast<ESymSolverStatus>(cuDSS_symbolic_factorization()) != SYMSOLVER_SUCCESS ) return SYMSOLVER_FATAL_ERROR;
        if (HaveIpData()) IpData().TimingStats().LinearSystemSymbolicFactorization().End();

        #if IPOPT_CUDSS_DEBUG == 1
        printf("cuDSS_symbolic_factorization END\n");
        #endif
    }

    initialized_ = true;

    return SYMSOLVER_SUCCESS;
}

Number *cuDSSSolverInterface::GetValuesArrayPtr()
{
    DBG_START_METH("cuDSSSolverInterface::GetValuesArrayPtr", dbg_verbosity);

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_get_matrix_values\n");
    #endif

    return cuDSS_get_matrix_values();
}

ESymSolverStatus cuDSSSolverInterface::MultiSolve(
    bool new_matrix, 
    const Index *ia, 
    const Index *ja, 
    Index nrhs, 
    Number *rhs_vals, 
    bool check_NegEVals, 
    Index numberOfNegEVals
)
{
    DBG_START_METH("cuDSSSolverInterface::MultiSolve", dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());

    if (new_matrix) {
        bool updated = cuDSS_update_matrix();
        ESymSolverStatus retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
        if ( retval != SYMSOLVER_SUCCESS ) return retval;
    }

    return Solve(ia, ja, nrhs, rhs_vals);
}

Index cuDSSSolverInterface::NumberOfNegEVals() const
{
    DBG_START_METH("cuDSSSolverInterface::NumberOfNegEVals", dbg_verbosity);
    DBG_ASSERT(negevals_ >= 0);
    
    return negevals_;
}

ESymSolverStatus cuDSSSolverInterface::Factorization(
    const Index *ia, 
    const Index *ja, 
    bool check_NegEVals, 
    Index numberOfNegEVals
)
{
    DBG_START_METH("cuDSSSolverInterface::Factorization", dbg_verbosity);

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_factorization START\n");
    #endif

    if (HaveIpData()) IpData().TimingStats().LinearSystemFactorization().Start();
    if ( static_cast<ESymSolverStatus>(cuDSS_factorization()) != SYMSOLVER_SUCCESS ) return SYMSOLVER_FATAL_ERROR;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "cuDSS Numerical Factorization.\n");
    if (HaveIpData()) IpData().TimingStats().LinearSystemFactorization().End();

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_factorization END\n");

    printf("cuDSS_get_inertia START\n");
    #endif

    negevals_ = Max(cuDSS_get_inertia(), numberOfNegEVals);

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_get_inertia END\n");
    #endif

    if( check_NegEVals && (numberOfNegEVals != negevals_) ) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Wrong inertia: required are %" IPOPT_INDEX_FORMAT ", but we got %" IPOPT_INDEX_FORMAT ".\n", numberOfNegEVals, negevals_);
        return SYMSOLVER_WRONG_INERTIA;
    }

    return SYMSOLVER_SUCCESS;
}

ESymSolverStatus cuDSSSolverInterface::Solve(
    const Index *ia, 
    const Index *ja, 
    Index nrhs, 
    Number *rhs_vals
)
{
    DBG_START_METH("cuDSSSolverInterface::Solve", dbg_verbosity);

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_solve START\n");
    #endif

    if (HaveIpData()) IpData().TimingStats().LinearSystemBackSolve().Start();
    if ( static_cast<ESymSolverStatus>(cuDSS_solve(nrhs, rhs_vals)) != SYMSOLVER_SUCCESS ) return SYMSOLVER_FATAL_ERROR;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "cuDSS Solve System.\n");
    if (HaveIpData()) IpData().TimingStats().LinearSystemBackSolve().End();

    #if IPOPT_CUDSS_DEBUG == 1
    printf("cuDSS_solve END\n");
    #endif

    return SYMSOLVER_SUCCESS;
}

} // namespace Ipopt