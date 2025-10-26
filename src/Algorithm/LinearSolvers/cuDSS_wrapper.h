#ifndef __CUDSS_WRAPPER_HPP__
#define __CUDSS_WRAPPER_HPP__

#ifdef CUDSS_SINGLE
typedef float Number;
#else
typedef double Number;
#endif

#ifdef CUDSS_INT64
typedef int64_t Index;
#else
typedef int Index;
#endif

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
int     cuDSS_get_inertia();

#endif