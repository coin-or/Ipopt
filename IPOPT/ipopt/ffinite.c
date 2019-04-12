/* Copyright (C) 2002, Carnegie Mellon University and others.
   All Rights Reserved.
   This code is published under the Common Public License.  */

/* $Id: ffinite.c 531 2004-03-11 01:31:07Z andreasw $  */
/*
    This routine implements the Fortran interface to the C function finite.
    Returns nonzero, if the argument is not Inf or NaN, returns nonzero
    if the argument is finite.
*/

/*
    Author:  Andreas Waechter     05/01/02  Release as Version IPOPT 2.0
*/

#ifdef __cplusplus
extern "C" {
#endif

#include <config.h>

/*  Make sure we do the correct casting from C int to Fortran INTEGER */
#include <Ipopt.h>

#ifdef HAVE_MATH_H
# include <math.h>
#endif
#ifdef HAVE_FLOAT_H
# include <float.h>
#endif
#ifdef HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

fint F77_FUNC(ffinite,FFINITE) (real *number)
{
  fint ret_val;

#ifdef MY_C_FINITE
  ret_val = (fint) MY_C_FINITE(*number);
#else
  ret_val = (fint) 1;
#endif

  return ret_val;
}

#ifdef __cplusplus
}
#endif
