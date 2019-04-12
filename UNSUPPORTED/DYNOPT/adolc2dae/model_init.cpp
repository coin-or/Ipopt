/*
  Copyright (C) 2002, Carnegie Mellon University and others.
  All Rights Reserved.
  This code is published under the Common Public License.
*/

/*  $Id: model_init.cpp 531 2004-03-11 01:31:07Z andreasw $  */
#include "adolc.h"
#include "f2c_adolc.h"

int model_(integer  *nz,
	   integer  *ny,
	   integer  *nu,
	   integer  *np,
	   adouble  *t,
	   adouble  *z,
	   adouble  *dz,
	   adouble  *y,
	   adouble  *u,
	   adouble  *p,
	   adouble  *f);

extern "C" {

  /* Call model once at the beginning to obtain tape */

#ifdef WIN32    // tell compiler that this function is to be exported
__declspec( dllexport )
#endif 
void model_init_
    (integer  *nz,
		integer  *ny,
		integer  *nu,
		integer  *np,
		double   *t,
		double   *z,
		double   *dz,
		double   *y,
		double   *u,
		double   *p,
		double   *f,
	  char     *flibname,
 	  int      libnamelen)
{ 
  int NZ = *nz;
  int NY = *ny;
  int NU = *nu;
  int NP = *np;


  /* Switch taping on */
  trace_on(1);

  adouble *tt, *zz, *dd, *yy, *uu, *pp, *ff;
  tt = new adouble[ 1 ];
  if( NZ > 0 )
    {
      zz = new adouble[NZ];
      dd = new adouble[NZ];
    }
  if( NY > 0 )
    yy = new adouble[NY];
  if( NU > 0 )
    uu = new adouble[NU];
  if( NP > 0 )
    pp = new adouble[NP];
  ff = new adouble[NZ+NY];

  int i;
  /* Select independent variables */
  *tt <<= *t;
  for( i=0; i < NZ; i++)
    zz[i] <<= z[i];
  for( i=0; i < NZ; i++)
    dd[i] <<= dz[i];
  for( i=0; i < NY; i++)
    yy[i] <<= y[i];
  for( i=0; i < NU; i++)
    uu[i] <<= u[i];
  for( i=0; i < NP; i++)
    pp[i] <<= p[i];

  /* Call model evaluation routine */
  model_(nz,ny,nu,np,tt,zz,dd,yy,uu,pp,ff);

  /* Select dependent variables */
  for( i=0; i < NZ+NY; i++)
    ff[i] >>= f[i];

  delete[] ff;
  if( NP > 0 )
    delete[] pp;
  if( NU > 0 )
    delete[] uu;
  if( NY > 0 )
    delete[] yy;
  if( NZ > 0 )
    {
      delete[] dd;
      delete[] zz;
    }
  delete[] tt;

  /*Switch off taping */
  trace_off();

}

}
