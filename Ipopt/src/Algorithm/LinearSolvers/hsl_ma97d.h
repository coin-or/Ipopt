/*
 * COPYRIGHT (c) 2012 The Science and Technology Facilities Council (STFC)
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Jonathan Hogg    STFC     2012-12-21
 */

#ifndef HSL_MA97D_H
#define HSL_MA97D_H

#ifndef ma97_default_control
#define ma97_control ma97_control_d
#define ma97_info ma97_info_d
#define ma97_default_control ma97_default_control_d
#define ma97_analyse ma97_analyse_d
#define ma97_analyse_coord ma97_analyse_coord_d
#define ma97_factor ma97_factor_d
#define ma97_factor_solve ma97_factor_solve_d
#define ma97_solve ma97_solve_d
#define ma97_free_akeep ma97_free_akeep_d
#define ma97_free_fkeep ma97_free_fkeep_d
#define ma97_finalise ma97_finalise_d
#define ma97_enquire_posdef ma97_enquire_posdef_d
#define ma97_enquire_indef ma97_enquire_indef_d
#define ma97_alter ma97_alter_d
#define ma97_solve_fredholm ma97_solve_fredholm_d
#define ma97_lmultiply ma97_lmultiply_d
#define ma97_sparse_fwd_solve ma97_sparse_fwd_solve_d
#endif

typedef double ma97pkgtype_d_;
typedef double ma97realtype_d_;

struct ma97_control_d {
    int f_arrays;             /* Use C or Fortran numbering */
    int action;               /* Continue on singularity if !=0 (true),
                                 otherwise abort */
    int nemin;                /* Supernode amalgamation if parent and child
                                 have fewer than nemin eliminations */
    ma97realtype_d_ multiplier;/* Amount of extra memory to allow for delays */
    int ordering;             /* Control scaling algorithm used:
                                 0 - user supplied order (order absent=identity)
                                 1 - AMD
                                 2 - MD (as in MA27)
                                 3 - METIS nested dissection
                                 4 - MA47
                                 5 - Automatic choice between 1 and 3 */
    int print_level;          /* <0 for no printing, 0 for basic, >1 for most */
    int scaling;              /* 0 user/none, 1 mc64, 2 mc77 */
    ma97realtype_d_ small;     /* Minimum value to count as non-zero */
    ma97realtype_d_ u;         /* Pivoting parameter */
    int unit_diagnostics;     /* Fortran unit for diagnostics (<0 disables) */
    int unit_error;           /* Fortran unit for error msgs (<0 disables) */
    int unit_warning;         /* Fortran unit for warning msgs (<0 disables) */
    long factor_min;          /* Min number of flops for parallel execution */
    int solve_blas3;          /* Use BLAS3 in solve in true, else BLAS2 */
    long solve_min;           /* Min number of entries for parallel exection */
    int solve_mf;             /* If true use m/f solve, else use s/n */
    ma97realtype_d_ consist_tol; /* Consistent equation tolerance */

    /* Reserve space for future interface changes */
    int ispare[5]; ma97realtype_d_ rspare[10];
};

struct ma97_info {
    int flag;                 /* <0 on error */
    int flag68;
    int flag77;
    int matrix_dup;           /* number duplicate entries in A */
    int matrix_rank;          /* matrix rank */
    int matrix_outrange;      /* number of out of range entries in A */
    int matrix_missing_diag;  /* number of zero diagonal entries in A */
    int maxdepth;             /* height of assembly tree */
    int maxfront;             /* maximum dimension of frontal matrix */
    int num_delay;            /* number of times a pivot was delayed */
    long num_factor;          /* number of entries in L */
    long num_flops;           /* number of floating point operations */
    int num_neg;              /* number of negative pivots */
    int num_sup;              /* number of supernodes in assembly tree */
    int num_two;              /* number of 2x2 pivots */
    int ordering;             /* ordering used (as per control.ordering) */
    int stat;                 /* error code from failed memory allocation */

    /* Reserve space for future interface changes */
    int ispare[5]; ma97realtype_d_ rspare[10];         
};

/* Set default values of control */
void ma97_default_control_d(struct ma97_control_d *control);
/* Perform symbolic analysis of matrix (sparse column entry) */
void ma97_analyse_d(int check, int n, const int ptr[], const int row[],
      ma97pkgtype_d_ val[], void **akeep, const struct ma97_control_d *control,
      struct ma97_info_d *info, int order[]);
/* Perform symbolic analysis of matrix (coordinate entry) */
void ma97_analyse_coord_d(int n, int ne, const int row[], const int col[],
      ma97pkgtype_d_ val[], void **akeep, const struct ma97_control_d *control,
      struct ma97_info_d *info, int order[]);
/* Perform numerical factorization, following call to ma97_analyse */
void ma97_factor_d(int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_d_ val[], void **akeep, void **fkeep,
      const struct ma97_control_d *control, struct ma97_info_d *info,
      ma97realtype_d_ scale[]);
/* Perform numerical factorization and solve, following call to ma97_analyse */
void ma97_factor_solve_d(int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_d_ val[], int nrhs, ma97pkgtype_d_ x[], int ldx,
      void **akeep, void **fkeep, const struct ma97_control_d *control,
      struct ma97_info_d *info, ma97realtype_d_ scale[]);
/* Perform forward and back substitutions, following call to ma97_factor */
void ma97_solve_d(int job, int nrhs, ma97pkgtype_d_ x[], int ldx,
      void **akeep, void **fkeep, const struct ma97_control_d *control,
      struct ma97_info_d *info);
/* Free memory in akeep */
void ma97_free_akeep_d(void **akeep);
/* Free memory in fkeep */
void ma97_free_fkeep_d(void **fkeep);
/* Free memory in akeep and fkeep */
void ma97_finalise_d(void **akeep, void **fkeep);
/* Return diagonal entries of L */
void ma97_enquire_posdef_d(void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      ma97realtype_d_ d[]);
/* Return diagonal, subdiagonal and/or pivot order of D */
void ma97_enquire_indef_d(void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      int *piv_order, ma97pkgtype_d_ *d);
/* Alter diagonal and subdiagonal of D */
void ma97_alter_d(const ma97pkgtype_d_ d[], void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info);
/* Fredholm alternative for singular systems */
void ma97_solve_fredholm_d(int nrhs,  int flag_out[], ma97pkgtype_d_ x[],
      int ldx, void **akeep, void **fkeep, const struct ma97_control_d *control,
      struct ma97_info_d *info);
/* Form (S^{-1}PL) X or (S^{-1}PL)^T X */
void ma97_lmultiply_d(int trans, int k, const ma97pkgtype_d_ x[], int ldx,
      ma97pkgtype_d_ y[], int ldy, void **akeep, void **fkeep,
      const struct ma97_control_d *control, struct ma97_info_d *info);
/* Perform a sparse forward solve */
void ma97_sparse_fwd_solve_d(int nbi, const int bindex[],
      const ma97pkgtype_d_ b[], const int order[], int *nxi, int xindex[],
      ma97pkgtype_d_ x[], void **akeep, void **fkeep,
      const struct ma97_control_d *control, struct ma97_info_d *info);

#endif
