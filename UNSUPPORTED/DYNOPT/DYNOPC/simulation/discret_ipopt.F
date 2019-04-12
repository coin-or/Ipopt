C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: discret_ipopt.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine DISCRET_IPOPT(TOL_DISCRET, NX, X, ZB, YB, UB, PB, IERR)
C
C     Solve the discretized model with fixed controls and parameters
C     with IPOPT.
C     For now: ignore bounds
C     Planned: give the option to move some/all of the controls and parameters
C              where deviation from initial point is penalized in objective
C              function
C     
C
C     Author:   Andreas Waechter    09-07-01
C
#ifdef WIN32
      use dflib
#endif

      implicit none
C
      double precision TOL_DISCRET ! TODO: make this the actual tolerance
                                   ! for IPOPT
      integer NX
      double precision X(NX)
      double precision ZB(2,*)
      double precision YB(2,*)
      double precision UB(2,*)
      double precision PB(2,*)
      integer IERR

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/

C     variables for IPOPT
      integer m
      integer nlb, nub
      integer ilb(NXMAX), iub(NXMAX)
      double precision bnds_l(NXMAX), bnds_u(NXMAX)
      double precision v_l(NXMAX), v_u(NXMAX), lam(NXMAX)
      integer iter

      double precision xbak(NXMAX)
      integer narg
      double precision args(50)
      character*20 cargs(50)

cCAW Right now do dynamic allocation
c      double precision rw(LRW_IPOPT)
c      integer          iw(LIW_IPOPT)
      integer PRW, PIW, MYALLOC

      integer CALC_NX

C     Store backup of DASPK simulation result in case IPOPT is not successful
      call DCOPY(NX, X, 1, xbak, 1)

C     No bounds for now
      nlb = 0
      nub = 0

C     For now: Do factorization based on first element
      CRIT_ELE = 1

C     For this use dummy (i.e. constant 0) objective function
      NZ_IN_OBJ = 0
      NY_IN_OBJ = 0
      NP_IN_OBJ = 0
      NE_IN_OBJ = 0

      call ADDCON_INIT(NZ, NY, 0, NCOL, 0, NAC)
!DEC$ ATTRIBUTES DLLIMPORT :: ADDCON_INIT
C     Set IPOPT parameters
      cargs(1) = 'ifull'
      args(1) = 0.d0
      cargs(2) = 'ihessvect'
      args(2) = 2.d0
      cargs(3) = 'dtol'
      args(3) = TOL_DISCRET
      cargs(4) = 'ierror'
      args(4) = -1.d0
      cargs(5) = 'itron2derivs'
      args(5) = 0.d0
      narg = 5
      call INITPARAMS(IERR, narg, args, cargs)
      if( IERR.ne.0 ) then
         write(*,90) IERR
 90      format(/,6x,'Error initializing IPOPT parameters. IERR =',i6)
         return
      endif

      write(*,100)
 100  format(/
     1     '>>>>> Calling IPOPT to fit data to discretization.'/)

C     Solve square problem with IPOPT
      m = CALC_NX(NZ, NY, 0, 0, NELE, NCOL)
CAW  Right now use dynamic allocation
c      call IPOPT(m, X, m, nlb, ilb, bnds_l, nub, iub, bnds_u, v_l, v_u,
c     1     lam, LRW_IPOPT, rw, LIW_IPOPT, iw, iter, IERR)
      PRW = MYALLOC(8*LRW_IPOPT)
      PIW = MYALLOC(4*LIW_IPOPT)
      call IPOPT(m, X, m, nlb, ilb, bnds_l, nub, iub, bnds_u, v_l, v_u,
     1     lam, LRW_IPOPT, %VAL(PRW), LIW_IPOPT, %VAL(PIW), iter, IERR)

      if( IERR.eq.0 ) then
         write(*,110)
 110     format(/
     1        '<<<<< IPOPT successfully obtained a feasible ',
     1        'starting point.'/)
      else
         write(*,120) IERR
 120     format(/
     1        6x,'IPOPT did not converge  --- IERR = ',i6,'.'/)
C     Restore simulation result
         call DCOPY(NX, xbak, 1, X, 1)
      endif

      return
      end
