/*
 * COPYRIGHT (c) 2011 Science and Technology Facilities Council (STFC)
 * Original date 20 September 2011
 * All rights reserved
 *
 * Written by: Jonathan Hogg
 *
 * Version 2.8.0
 *
 * THIS FILE ONLY may be redistributed under the below modified BSD licence.
 * All other files distributed as part of the HSL_MA97 package
 * require a licence to be obtained from STFC and may NOT be redistributed
 * without permission. Please refer to your licence for HSL_MA97 for full terms
 * and conditions. STFC may be contacted via hsl(at)stfc.ac.uk.
 *
 * Modified BSD licence (this header file only):
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of STFC nor the names of its contributors may be used
 *    to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL STFC BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef HSL_MA97S_H
#define HSL_MA97S_H

#ifndef ma97_default_control
#define ma97_control ma97_control_s
#define ma97_info ma97_info_s
#define ma97_default_control ma97_default_control_s
#define ma97_analyse ma97_analyse_s
#define ma97_analyse_coord ma97_analyse_coord_s
#define ma97_factor ma97_factor_s
#define ma97_factor_solve ma97_factor_solve_s
#define ma97_solve ma97_solve_s
#define ma97_free_akeep ma97_free_akeep_s
#define ma97_free_fkeep ma97_free_fkeep_s
#define ma97_finalise ma97_finalise_s
#define ma97_enquire_posdef ma97_enquire_posdef_s
#define ma97_enquire_indef ma97_enquire_indef_s
#define ma97_alter ma97_alter_s
#define ma97_solve_fredholm ma97_solve_fredholm_s
#define ma97_lmultiply ma97_lmultiply_s
#define ma97_sparse_fwd_solve ma97_sparse_fwd_solve_s
#endif

typedef float ma97pkgtype_s_;
typedef float ma97realtype_s_;

struct ma97_control_s {
    int f_arrays;             /* Use C or Fortran numbering */
    int action;               /* Continue on singularity if !=0 (true),
                                 otherwise abort */
    int nemin;                /* Supernode amalgamation if parent and child
                                 have fewer than nemin eliminations */
    ma97realtype_s_ multiplier;/* Amount of extra memory to allow for delays */
    int ordering;             /* Control scaling algorithm used:
                                 0 - user supplied order (order absent=identity)
                                 1 - AMD
                                 2 - MD (as in MA27)
                                 3 - METIS nested dissection
                                 4 - MA47
                                 5 - Automatic choice between 1 and 3 */
    int print_level;          /* <0 for no printing, 0 for basic, >1 for most */
    int scaling;              /* 0 user/none, 1 mc64, 2 mc77 */
    ma97realtype_s_ small;     /* Minimum value to count as non-zero */
    ma97realtype_s_ u;         /* Pivoting parameter */
    int unit_diagnostics;     /* Fortran unit for diagnostics (<0 disables) */
    int unit_error;           /* Fortran unit for error msgs (<0 disables) */
    int unit_warning;         /* Fortran unit for warning msgs (<0 disables) */
    long factor_min;          /* Min number of flops for parallel execution */
    int solve_blas3;          /* Use BLAS3 in solve in true, else BLAS2 */
    long solve_min;           /* Min number of entries for parallel exection */
    int solve_mf;             /* If true use m/f solve, else use s/n */
    ma97realtype_s_ consist_tol; /* Consistent equation tolerance */

    /* Reserve space for future interface changes */
    int ispare[5]; ma97realtype_s_ rspare[10];
};

struct ma97_info_s {
    int flag;                 /* <0 on error */
    int flag68;
    int flag77;
    int matrix_dup;           /* number duplicate entries in A */
    int matrix_rank;          /* matrix rank */
    int matrix_outrange;      /* number of out of range entries in A */
    int matrix_missing_diag;  /* number of zero diagonal entries in A */
    int maxdepth;             /* height of assembly tree */
    int maxfront;             /* maximum no. rows in a supernode */
    int num_delay;            /* number of times a pivot was delayed */
    long num_factor;          /* number of entries in L */
    long num_flops;           /* number of floating point operations */
    int num_neg;              /* number of negative pivots */
    int num_sup;              /* number of supernodes in assembly tree */
    int num_two;              /* number of 2x2 pivots */
    int ordering;             /* ordering used (as per control.ordering) */
    int stat;                 /* error code from failed memory allocation */
    int maxsupernode;         /* maximum no. columns in a supernode */

    /* Reserve space for future interface changes */
    int ispare[4]; ma97realtype_s_ rspare[10];
};

/* Set default values of control */
void ma97_default_control_s(struct ma97_control_s *control);
/* Perform symbolic analysis of matrix (sparse column entry) */
void ma97_analyse_s(int check, int n, const int ptr[], const int row[],
      ma97pkgtype_s_ val[], void **akeep, const struct ma97_control_s *control,
      struct ma97_info_s *info, int order[]);
/* Perform symbolic analysis of matrix (coordinate entry) */
void ma97_analyse_coord_s(int n, int ne, const int row[], const int col[],
      ma97pkgtype_s_ val[], void **akeep, const struct ma97_control_s *control,
      struct ma97_info_s *info, int order[]);
/* Perform numerical factorization, following call to ma97_analyse */
void ma97_factor_s(int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_s_ val[], void **akeep, void **fkeep,
      const struct ma97_control_s *control, struct ma97_info_s *info,
      ma97realtype_s_ scale[]);
/* Perform numerical factorization and solve, following call to ma97_analyse */
void ma97_factor_solve_s(int matrix_type, const int ptr[], const int row[],
      const ma97pkgtype_s_ val[], int nrhs, ma97pkgtype_s_ x[], int ldx,
      void **akeep, void **fkeep, const struct ma97_control_s *control,
      struct ma97_info_s *info, ma97realtype_s_ scale[]);
/* Perform forward and back substitutions, following call to ma97_factor */
void ma97_solve_s(int job, int nrhs, ma97pkgtype_s_ x[], int ldx,
      void **akeep, void **fkeep, const struct ma97_control_s *control,
      struct ma97_info_s *info);
/* Free memory in akeep */
void ma97_free_akeep_s(void **akeep);
/* Free memory in fkeep */
void ma97_free_fkeep_s(void **fkeep);
/* Free memory in akeep and fkeep */
void ma97_finalise_s(void **akeep, void **fkeep);
/* Return diagonal entries of L */
void ma97_enquire_posdef_s(void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      ma97realtype_s_ d[]);
/* Return diagonal, subdiagonal and/or pivot order of D */
void ma97_enquire_indef_s(void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info,
      int *piv_order, ma97pkgtype_s_ *d);
/* Alter diagonal and subdiagonal of D */
void ma97_alter_s(const ma97pkgtype_s_ d[], void **akeep, void **fkeep,
      const struct ma97_control *control, struct ma97_info *info);
/* Fredholm alternative for singular systems */
void ma97_solve_fredholm_s(int nrhs, int flag_out[], ma97pkgtype_s_ x[],
      int ldx, void **akeep, void **fkeep, const struct ma97_control_s *control,
      struct ma97_info_s *info);
/* Form (S^{-1}PL) X or (S^{-1}PL)^T X */
void ma97_lmultiply_s(int trans, int k, const ma97pkgtype_s_ x[], int ldx,
      ma97pkgtype_s_ y[], int ldy, void **akeep, void **fkeep,
      const struct ma97_control_s *control, struct ma97_info_s *info);
/* Perform a sparse forward solve */
void ma97_sparse_fwd_solve_s(int nbi, const int bindex[],
      const ma97pkgtype_s_ b[], const int order[], int *nxi, int xindex[],
      ma97pkgtype_s_ x[], void **akeep, void **fkeep,
      const struct ma97_control_s *control, struct ma97_info_s *info);

#endif
