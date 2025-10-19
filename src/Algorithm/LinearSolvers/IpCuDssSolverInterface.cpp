// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpCuDssSolverInterface.hpp"

void Ipopt::cuDSSSolverInterface::RegisterOptions(
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
}