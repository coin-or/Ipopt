/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/* f2c.h  --  This header file has been modified in order to work with ADOL-C
   03-31-00   Andreas Waechter, Dept. Chem. Eng., Carnegie Mellon University
   $Id: f2c_adolc.h 531 2004-03-11 01:31:07Z andreasw $  */

#ifndef __F2C_ADOLC__
#define __F2C_ADOLC__
#include "adolc.h"

#ifdef __osf__
typedef int integer;
#else
typedef long int integer;
#endif
typedef adouble doublereal;
typedef long int logical;

/* AW: My definitions: */
#define pow_dd(x,y) (pow(*x,*y))
#define dmax1(a,b) (fmax(a,b))
#define dmin1(a,b) (fmin(a,b))
#define condassign_(a,b,c,d) (condassign(*a,*b,*c,*d))
#define log(a) log( (adouble) a)

#define abs(x) ((x) >= 0 ? (x) : -(x))
#ifdef min
#undef min
#endif
#define min(a,b) ((a) <= (b) ? (a) : (b))
#ifdef max
#undef max
#endif
#define max(a,b) ((a) >= (b) ? (a) : (b))
/*
#define dabs(x) (doublereal)abs(x)
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
*/
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))

#endif

