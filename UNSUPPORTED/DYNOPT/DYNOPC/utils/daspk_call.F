C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: daspk_call.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine DASPK_CALL(TASK, TOL_DASPK, ZB, Z0, DZ, YB, Y0, UB,
     1     U0, PB, P0, X, XIN, NX, CNR_DASPK, IERR)
C
C     Interface to DASPK; performs simulation of process using DASPK between
C     the finite elements
C
C     Authors:  Yi-dong Lang
C               Andreas Waechter    09-06-01  rewritten / replaced DASSL
C
      implicit none
C
      integer TASK              ! =0: Called in simulation - set controls
                                !     constant
                                ! =1: Called after optimization - values
                                !     of control variables according to optimal
                                !     solution
      double precision TOL_DASPK
      double precision ZB(2,*), Z0(*), DZ(*)
      double precision YB(2,*), Y0(*)
      double precision UB(2,*), U0(*)
      double precision PB(2,*), P0(*)
      double precision X(*)
      double precision XIN(*)   ! only for TASK = 1: optimal X vector
      integer NX                ! Size of discretization (output!)
      integer CNR_DASPK         ! Fortran Channel number for DASPK 
      integer IERR

      include 'DYNAUX.INC'
      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DAENLP/

C     Variables for DASPK
      integer neq
      double precision t
      double precision yd(NZMAX+NYMAX)
      double precision yprime(NZMAX+NYMAX)
      double precision tout
      integer info(20)
      double precision rtol, atol
      integer idid
      integer lrw
      parameter( lrw = 50 + 9*(NZMAX+NYMAX) + (NZMAX+NYMAX)**2 )
      double precision rwork(lrw)
      integer liw
      parameter( liw = 40 + 2*(NZMAX+NYMAX) )
      integer iwork(liw)
      external RES_DASPK, JAC_DASPK

      integer i, lele, lcol
      double precision dti
C
C     First we need to set all control variables that do not have a specified
C     profile as fixed for the DAEMODEL interface - same for parameters
C     (see DYNAUX.INC)
C
      if( TASK.eq.0 )then
         if( NU_PROF.ne.NU ) then
            do i = 1, NU
               if( IU_PROF(i).eq.0 ) then
                  IU_PROF(i) = 1
                  call DCOPY(NELE+1, U0(i), 0, U_PROF(i,1), NUMAX)
               endif
            enddo
            NU_PROF = NU
         endif
         if( NP_FIX.ne.NP ) then
            do i = 1, NP
               if( IP_FIX(i).eq.0 ) then
                  IP_FIX(i) = 1
                  P_FIX(i) = P0(i)
               endif
            enddo
            NP_FIX = NP
         endif
      else
         do i = 1, NU
            IU_PROF(i) = 1
         enddo
         NU_PROF = NU
         do i = 1, NP
            IP_FIX(i) = 1
         enddo
         NP_FIX = NP
      endif
C
C     Prepare parameters for DASPK
C
      neq = NZ + NY
      NX = 0

C     Use default for most INFO parameters
      do i = 1, 20
         info(i) = 0
      enddo
      info( 5) = 1              ! use Jacobian routine
      if( CNR_DASPK.gt.0 ) then
         info(18) = 2           ! for now, we want some debug info from DASPK
         CALL XSETUN(CNR_DASPK) ! Set output unit number for DASPK
      endif

      rtol = 0.d0
      atol = TOL_DASPK
C ============================================================================
C
C     Obtain consistent initial conditions from DASPK
C
C ============================================================================


C     Prepare initial point for DASPK
      call DCOPY(NZ, Z0  , 1, yd    (   1), 1)
      call DCOPY(NY, Y0  , 1, yd    (NZ+1), 1)
      call DCOPY(NZ, DZ  , 1, yprime(   1), 1)
      call DCOPY(NY, 0.d0, 0, yprime(NZ+1), 1)

      write(*,90)
 90   format(/
     1     '>>>>> Starting simulation of process with DASPK .'/)

      t= TI(1)

      if( TASK.eq.0 ) goto 120  ! perform initialization with DASPK only
                                ! for TASK = 1

C     We want DASPK to make initial conditions consistent:
      info(11) = 1
      do i = 40+1, 40+NZ
         iwork(i) = 1           ! first NZ variables in yd are differentiable
      enddo
      do i = 40+NZ+1, 40+NZ+NY
         iwork(i) = -1          ! remaining NY variables in yb are algebraic
      enddo
      info(14) = 1              ! we want consistent initial conditions back
                                ! from DASPK before starting the simulation

      tout = TI(2)
      call DDASPK(RES_DASPK, neq, t, yd, yprime, tout, info, rtol, atol,
     1     idid, rwork, lrw, iwork, liw, XIN, TASK, JAC_DASPK,
     1     JAC_DASPK)

      if( idid.eq.4 ) then
         write(*,100)
 100     format(/
     1     6x,'DASPK successfully determined ',
     2        'consistent initial conditions'/)
      else
         write(*,101) idid
 101     format(/
     1     6x,'***********************************************'/
     2     6x,'*     Problem with during initialization!     *'/
     3     6x,'*         DASPK returned IDID = ',i4,'          *'/
     4     6x,'* Still need to decide what to do.  Now abort *'/
     4     6x,'***********************************************'/)
         IERR = 1
         return
      endif

 120  continue
C     Copy initial conditions into X
      call DCOPY(NZ, yd, 1, X(NX+1), 1)
      NX = NX + NZ

C ============================================================================
C
C     Simulate the process with DASPK
C
C ============================================================================

C     Make sure it is not again trying to determine the initial conditions
      info (11) = 0

C     call the initialization routine for DAE2NLP, in order to have RHO
C     in DAE2NLP.INC right
      call INITD2N(NCOL)

      do lele = 1, NELE
         dti = TI(lele+1) - TI(lele)
         do lcol = 1, NCOL
C
            tout = TI(lele) + dti*RHO(lcol)
            call DDASPK(RES_DASPK, neq, t, yd, yprime, tout, info,
     1           rtol, atol, idid, rwork, lrw, iwork, liw, XIN, TASK,
     2           JAC_DASPK, JAC_DASPK)
#ifndef WIN32
            call CPRINTF('%s','*')
#else
            write(*,1001) '*'
 1001       format(a1\)
#endif
            if( idid.ne.3 ) then ! something is wrong
               write(*,*) 'Error in DASPK; IDID = ',idid
               IERR = idid
               return
            endif
C
C     At this point later bound check
C

C     Copy result into X
            call DCOPY(NZ, yprime, 1, X(NX+1), 1)
            NX = NX + NZ
            call DCOPY(NY, yd(NZ+1), 1, X(NX+1), 1)
            NX = NX + NY
C
         enddo

C     If for some reason (depending on collocation points) T <> TI(i+1)
C     do the remaining step with DDASPK
         if( tout.lt.TI(lele+1) ) then
            tout = TI(lele+1)
            call DDASPK(RES_DASPK, neq, t, yd, yprime, tout, info,
     1           rtol, atol, idid, rwork, lrw, iwork, liw, XIN, TASK,
     2           JAC_DASPK, JAC_DASPK)
            if( idid.ne.3 ) then ! something is wrong
               write(*,*) 'Error in DASPK; IDID = ',idid
               IERR = idid
               return
            endif
         endif

C     Copy Z into X
         call DCOPY(NZ, yd, 1, X(NX+1), 1)
         NX = NX + NZ
      enddo

C ============================================================================
C
C     Write profiles of controls and parameters into X
C
C ============================================================================

      do lele = 1, NELE
         dti = TI(lele+1) - TI(lele)
         do lcol = 1, NCOL
            call DCOPY(NU, U_PROF(1,lele), 1, X(NX+1), 1)
            call DSCAL(NU, (1.d0-RHO(lcol)), X(NX+1), 1)
            call DAXPY(NU, RHO(lcol), U_PROF(1,lele+1), 1, X(NX+1), 1)
            NX = NX + NU
         enddo
      enddo

      call DCOPY(NP, P_FIX, 1, X(NX+1), 1)
      NX = NX + NP

      write(*,110)
 110  format(//
     1     '<<<<< DASPK successfully simulated the process.'/)

      return
      end

C ============================================================================

C $Id: daspk_call.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine RES_DASPK(T, YD, YPRIME, CJ, DELTA, IRES, XIN, TASK)
C
C     Compute residual for DASPK using DAEMODEL_F.
C
C     Author:     Andreas Waechter     09-06-01
C
      implicit none
      double precision T
      double precision YD(*)
      double precision YPRIME(*)
      double precision CJ
      double precision DELTA(*)
      integer IRES
      double precision XIN(*)
      integer TASK

C
      include 'DYNAUX.INC'
C
      double precision dummy(1)
      integer iele, GET_IELE4T
      double precision zval(NZMAX), dzval(NZMAX), yval(NYMAX)

      if( TASK.eq.1 ) then
C
C     Need to make sure that the correct value for control is used!
C
         iele = GET_IELE4T(T)
         call APPSLN(T, NZ, NY, NU, NELE, NCOL, TI, XIN,
     1     zval, dzval, yval, U_PROF(1,iele))
         if( iele.le.NELE ) then
            call DCOPY(NU, U_PROF(1,iele), 1, U_PROF(1,iele+1), 1)
         else
            call DCOPY(NU, U_PROF(1,iele), 1, U_PROF(1,iele-1), 1)
         endif
      endif

      call DAEMODEL_F( NZ, NY, 0, 0, T, YD(1),
     1     YPRIME(1), YD(NZ+1), dummy, dummy, DELTA )
#ifndef WIN32
      call CPRINTF('%s','.')
#else
      write(*,100) '.'
  100 format(a1\)
#endif
      return
      end

C ============================================================================

C $Id: daspk_call.F 531 2004-03-11 01:31:07Z andreasw $
      subroutine JAC_DASPK(T, YD, YPRIME, PD, CJ, XIN, TASK)
C
C     Compute Jacobian for DASPK using DAEMODEL_DF.
C
C     Author:     Andreas Waechter     09-06-01
C
      implicit none
      double precision T
      double precision YD(*)
      double precision YPRIME(*)
      double precision CJ
      double precision PD(*)
      double precision XIN(*)
      integer TASK
C
      include 'DYNAUX.INC'
C
      double precision dummy(1)
      double precision buffer((NZMAX+NYMAX)*(2*NZMAX+NYMAX))
                                ! buffer for Jacobian from DAEMODEL_DF
                                ! in JAC_DASPK
      integer GET_IELE4T, iele
      double precision zval(NZMAX), dzval(NZMAX), yval(NYMAX)

      if( TASK.eq.1 ) then
C
C     Need to make sure that the correct value for control is used!
C
         iele = GET_IELE4T(T)
         call APPSLN(T, NZ, NY, NU, NELE, NCOL, TI, XIN,
     1     zval, dzval, yval, U_PROF(1,iele))
         if( iele.le.NELE ) then
            call DCOPY(NU, U_PROF(1,iele), 1, U_PROF(1,iele+1), 1)
         else
            call DCOPY(NU, U_PROF(1,iele), 1, U_PROF(1,iele-1), 1)
         endif
      endif

C     Get Jacobian RPAR in DAEMODEL_DF order
      call DAEMODEL_DF( NZ, NY, 0, 0, T, YD(1), YPRIME(1),
     1     YD(nz+1), dummy, dummy, buffer, NZ+NY )

C     Construct Jacobian PD as DASPK wants it
      call DCOPY((NZ+NY)*NZ, buffer(1), 1, PD(1), 1)
      call DCOPY((NZ+NY)*NY, buffer((NZ+NY)*2*NZ+1), 1,
     1     PD((NZ+NY)*NZ+1), 1)
      call DAXPY((NZ+NY)*NZ, CJ, buffer((NZ+NY)*NZ+1), 1, PD(1), 1)

      return
      end
