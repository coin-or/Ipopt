/*
 * COPYRIGHT (c) 2011 The Science and Technology Facilities Council (STFC)
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Jonathan Hogg    STFC     2011-05-18
 */

#ifndef HSL_MA77D_H
#define HSL_MA77D_H

#ifndef ma77_default_control
#define ma77_control ma77_control_d
#define ma77_info ma77_info_d
#define ma77_default_control ma77_default_control_d
#define ma77_open_nelt ma77_open_nelt_d
#define ma77_open ma77_open_d
#define ma77_input_vars ma77_input_vars_d
#define ma77_input_reals ma77_input_reals_d
#define ma77_analyse ma77_analyse_d
#define ma77_factor ma77_factor_d
#define ma77_factor_solve ma77_factor_solve_d
#define ma77_solve ma77_solve_d
#define ma77_resid ma77_resid_d
#define ma77_scale ma77_scale_d
#define ma77_enquire_posdef ma77_enquire_posdef_d
#define ma77_enquire_indef ma77_enquire_indef_d
#define ma77_alter ma77_alter_d
#define ma77_restart ma77_restart_d
#define ma77_finalise ma77_finalise_d
#endif

typedef double ma77pkgtype_d_;

/* Data type for user controls */
struct ma77_control_d {
   /* Note: 0 is false, non-zero is true */

   /* C/Fortran interface related controls */
   int f_arrays; /* Treat arrays as 1-based (Fortran) if true or 0-based (C) if
                    false. */

   /* Printing controls */
   int print_level;
   int unit_diagnostics;   /* unit for diagnostic messages
                              Printing is suppressed if unit_diagnostics < 0. */
   int unit_error;         /* unit for error messages
                              Printing is suppressed if unit_error  <  0.     */
   int unit_warning;       /* unit for warning messages
                              Printing is suppressed if unit_warning  <  0.   */

   /* Controls used by MA77_open */
   int bits;
   int buffer_lpage[2];
   int buffer_npage[2];
   long int file_size;
   long int maxstore;
   long int storage[3];

   /* Controls used by MA77_analyse */
   int nemin;  /* Node amalgamation parameter. A child node is merged with its
                  parent if they both involve fewer than nemin eliminations.*/

   /* Controls used by MA77_scale */
   int maxit;
   int infnorm;
   ma77pkgtype_d_ thresh;

   /* Controls used by MA77_factor with posdef true */
   int nb54;

   /* Controls used by MA77_factor with posdef false */
   int action;    /* Keep going even if matrix is singular if true, or abort
                     if false */
   ma77pkgtype_d_ multiplier;
   int nb64;
   int nbi;
   ma77pkgtype_d_ small;
   ma77pkgtype_d_ static_;
   long int storage_indef;
   ma77pkgtype_d_ u;       /* Pivot tolerance*/
   ma77pkgtype_d_ umin;    /* Minimum pivot tolerance*/
};

/***************************************************/

/* data type for returning information to user.*/
struct ma77_info_d {
   ma77pkgtype_d_ detlog;
   int detsign;
   int flag;
   int iostat;
   int matrix_dup;
   int matrix_rank;
   int matrix_outrange;
   int maxdepth;
   int maxfront;
   long int minstore;
   int ndelay;
   long int nfactor;
   long int nflops;
   int niter;
   int nsup;
   int num_neg;
   int num_nothresh;
   int num_perturbed;
   int ntwo;
   int stat;
   int index[4];
   long int nio_read[2];
   long int nio_write[2];
   long int nwd_read[2];
   long int nwd_write[2];
   int num_file[4];
   long int storage[4];
   int tree_nodes;
   int unit_restart;
   int unused;
   ma77pkgtype_d_ usmall;
};

/* Initialise control with default values */
void ma77_default_control(struct ma77_control_d *control);
void ma77_open_nelt(const int n, const char* fname1, const char* fname2,
   const char *fname3, const char *fname4, void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const int nelt);
void ma77_open(const int n, const char* fname1, const char* fname2,
   const char *fname3, const char *fname4, void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info);
void ma77_input_vars(const int idx, const int nvar, const int list[],
   void **keep, const struct ma77_control_d *control, struct ma77_info_d *info);
void ma77_input_reals(const int idx, const int length,
   const ma77pkgtype_d_ reals[], void **keep, const struct ma77_control_d *control,
   struct ma77_info_d *info);
/* Analyse the sparsity pattern and prepare for factorization */
void ma77_analyse(const int order[], void **keep,
   const struct ma77_control_d *control, struct ma77_info_d *info);
/* To factorize the matrix */
void ma77_factor(const int posdef, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale);
/* To factorize the matrix AND solve AX = B */
void ma77_factor_solve(const int posdef, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale, const int nrhs, const int lx,
   ma77pkgtype_d_ rhs[]);
/* To solve AX = B using the computed factors */
void ma77_solve(const int job, const int nrhs, const int lx, ma77pkgtype_d_ x[],
   void **keep, const struct ma77_control_d *control, struct ma77_info_d *info,
   const ma77pkgtype_d_ *scale);
void ma77_resid(const int nrhs, const int lx, const ma77pkgtype_d_ x[],
   const int lresid, ma77pkgtype_d_ resid[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   ma77pkgtype_d_ *anorm_bnd);
void ma77_scale(ma77pkgtype_d_ scale[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info,
   ma77pkgtype_d_ *anorm);
void ma77_enquire_posdef(ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info);
void ma77_enquire_indef(int piv_order[], ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info);
void ma77_alter(const ma77pkgtype_d_ d[], void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info);
void ma77_restart(const char *restart_file, const char *fname1, 
   const char *fname2, const char *fname3, const char *fname4, void **keep, 
   const struct ma77_control_d *control, struct ma77_info_d *info);
/* To clean up memory in keep */
void ma77_finalise(void **keep, const struct ma77_control_d *control,
   struct ma77_info_d *info);

#endif
