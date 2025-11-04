// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Antonio Cioffi                          2025-10-19
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __CUDSS_WRAPPER_HPP__
#define __CUDSS_WRAPPER_HPP__

#ifdef CUDSS_SINGLE
typedef float Number;
#else
typedef double Number;
#endif

typedef int64_t Index;

struct cuDSS_config_settings
{
    int algReorder;
    int algFactor;
    int algPivotEps;
    int useMatching;
    int algMatching;
    int nIterSteps;
    int pivotType;
    double pivotThr;
    double pivotEps;
    Index maxLUnnz;
    int nThreads;
    int ndNLevels;
    int useSP;
    int schurMode;
    int deterministic;
};

void    cuDSS_initialize();
void    cuDSS_terminate();
bool    cuDSS_config_create_and_set(cuDSS_config_settings);
void    cuDSS_initialize_structure(Index, Index, const Index*, const Index*);
int     cuDSS_reordering();
int     cuDSS_symbolic_factorization();
Number* cuDSS_get_matrix_values();
int     cuDSS_factorization();
int     cuDSS_refactorization();
int     cuDSS_solve(Index, Number*);
Index   cuDSS_get_inertia();
bool    cuDSS_update_matrix();

#endif