/*
 * COPYRIGHT (c) 2011 Science and Technology Facilities Council (STFC)
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Jonathan Hogg    STFC     2011-02-25
 */

#ifndef HSL_MC68I_H
#define HSL_MC68I_H

struct mc68_control {
   /* Extra options for C version */
   int f_array_in;      /* 0 for C array indexing, 1 for Fortran indexing */
   int f_array_out;     /* 0 for C array indexing, 1 for Fortran indexing
                         * NOTE: 2x2 pivot information discarded if C indexing
                         * is used for output! */
   int min_l_workspace; /* Initial size of workspace, as argument in Fortran */
   /* Options from Fortran version */
   int lp;              /* stream number for error messages */
   int wp;              /* stream number for warning messages */
   int mp;              /* stream number for diagnostic messages */
   int nemin;           /* stream number for diagnostic messages */
   int print_level;     /* amount of informational output required */
   int row_full_thresh; /* percentage threshold for full row */
   int row_search;      /* Number of rows searched for pivot with ord=6 */
};

struct mc68_info {
   int flag;            /* error/warning flag */
   int iostat;          /* holds Fortran iostat parameter */
   int stat;            /* holds Fortran stat parameter */
   int out_range;       /* holds number of out of range entries ignored */
   int duplicate;       /* holds number of duplicate entries */
   int n_compressions;  /* holds number of compressions in order */
   int n_zero_eigs;     /* holds the number of zero eigs from ma47 */
   int l_workspace;     /* holds length of workspace iw used in order */
   int zb01_info;       /* holds flag from zb01_expand1 call */
   int n_dense_rows;    /* holds number of dense rows from amdd */
};

/* Set default values for control struct */
void mc68_default_control(struct mc68_control *control);
/* Perform ordering */
void mc68_order(const int ord, const int n, const int ptr[], const int row[],
   int perm[], const struct mc68_control *control, struct mc68_info *info);

#endif
