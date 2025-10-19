// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpCuDssSolverInterface.hpp"

namespace Ipopt
{

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
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_ubatch_size",
        "The number of matrices in a uniform batch of systems to be processed by cuDSS.",
        1,
        1);
    roptions->AddLowerBoundedIntegerOption(
        "cuDSS_ubatch_index",
        "-1 or a 0-based index of matrix in a uniform batch which will be processed during factorization or solve phase.",
        -1,
        -1);
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
    return true;
}

} // namespace Ipopt