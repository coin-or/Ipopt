C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: DYNAUX.INC 531 2004-03-11 01:31:07Z andreasw $
C     Parameters and common block for DYNOPT/auxil
C
C     Author:  Andreas Waechter     09-03-01
C

C     Size of work space for IPOPT
      integer LIW_IPOPT, LRW_IPOPT
CHUGE      parameter( LIW_IPOPT = 130000000, LRW_IPOPT = 200000000 )
CMED      parameter( LIW_IPOPT = 20000000, LRW_IPOPT = 20000000 )
CSMALL      parameter( LIW_IPOPT = 2000000, LRW_IPOPT = 2000000 )
      parameter( LIW_IPOPT = 2000000, LRW_IPOPT = 2000000 )
C
C     maximal size of Z, Y, U, P  vectors
C
      integer NZMAX, NYMAX, NUMAX, NPMAX
CHUGE      parameter( NZMAX = 600, NYMAX = 3000, NUMAX = 20, NPMAX = 20)
      parameter( NZMAX = 60, NYMAX = 300, NUMAX = 10, NPMAX = 10)
C
C     maximum number of finite elements and collocation points
C
      integer NCOLMAX, NELEMAX
CHUGE      parameter( NCOLMAX = 6, NELEMAX = 2000 )
      parameter( NCOLMAX = 6, NELEMAX = 200 )
C
C     maximal size of vector X
C
      integer NXMAX
      parameter( NXMAX = 2500000 )
C
C     "very large number", i.e. if bounds are beyond -VLARGE and VLARGE,
C     respectively, they are considered to be infinitiy.
C
      double precision VLARGE
      parameter( VLARGE = 1.d20 )
C
C     TI:    time discretization
C     ZINIT: initial conditions for differential variables
C
      double precision TI(NELEMAX+1), ZINIT(NZMAX)
C
C     NZ:       number of differential variables
C     NY:       number of algebraic variables
C     NU:       number of control variables
C     NP:       number of parameters
C
C     NOTE:  NZ, NY, NU, NP might differ from what the optimizer sees in case
C            there are given profiles for U's or fixed values for R's (see
C            below)
C
C     NELE:     number of finite elements
C     NCOL:     number of collocation points
C     NAC:      number of additional constraints (see dae2nlp.f)
C     CRIT_ELE: number of element to be used for basis partitioning
C                    (see dae2nlp.f)
C
      integer NZ, NY, NU, NP, NELE, NCOL, NAC, CRIT_ELE
C
C     Information about the control variables with given profile.
C     The way implemented now is not very memory efficient...
C
      integer NU_PROF           ! number of controls with given profile
      integer IU_PROF(NUMAX)    ! for each control i, IU_PROF(i) is non-zero
                                ! if its profile is given in U_PROF(i,...)
      double precision U_PROF(NUMAX, NELEMAX+1)
                                ! values of those controls at element
                                ! boundaries (might take a lot of memory...)
C
C     Information about fixed parameters
C
      integer NP_FIX            ! number of fixed parameters
      integer IP_FIX(NPMAX)     ! for each parameter i, IP_FIX(i) is non-zero
                                ! if its value is given in P_FIX(i)
      double precision P_FIX(NPMAX)
                                ! values of those parameters

C     Counter for model_f and model_df  calls
      integer NO_MODEL_F_CALL, NO_MODEL_DF_CALL


      common /DYNAUX/ TI, ZINIT, U_PROF, P_FIX, IU_PROF, IP_FIX,
     1     NU_PROF, NP_FIX, NZ, NY, NU, NP, NELE, NCOL, NAC, CRIT_ELE,
     2     NO_MODEL_F_CALL, NO_MODEL_DF_CALL
      save   /DYNAUX/
