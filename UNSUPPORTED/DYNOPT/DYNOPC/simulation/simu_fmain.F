C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C =============================================================================
C $Id: simu_fmain.F 531 2004-03-11 01:31:07Z andreasw $
C =============================================================================
      subroutine FMAIN( MODNAM, STUBNAM )
C =============================================================================
C
C     Main program for the Simulation
C
C     Authors:  Yi-dong Lang
C               Andreas Waechter    09-04-01   cleaned up / made portable
C
C =============================================================================

#ifdef WIN32
      use dflib
#endif
      implicit none
      character*(*) MODNAM
      character*(*) STUBNAM

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
#ifdef WIN32
      include 'DYNGRA.INC'
      external PRINTSUB
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
#endif

C     Bounds and initial values for model variables read from .cmp file
      double precision zb(2,NZMAX), z0(NZMAX)
      double precision yb(2,NYMAX), y0(NYMAX)
      double precision ub(2,NUMAX), u0(NUMAX)
      double precision pb(2,NPMAX), p0(NPMAX)

C     Derivatives for differential variables
      double precision dz(NZMAX)

C     Dummy residual for initialization of ADOL-C
      double precision f(NZMAX+NYMAX)

C     Tolerances read from .cmp file
      double precision tol_daspk, tol_discret, tol_opt

C     Vector X for discretization
      double precision x(NXMAX)
      integer nx                ! actual size of X

C     Specification of model executable
      character*512 filename

      integer i, cnr_daspk, cnr_nleq1s, ierr, nlines
      double precision tol_nleq1s, dummy
      character*80 lines(7)
      logical ex

#ifdef WIN32
      integer(4) event, ix, iy, i4, i4k
      type(qwinfo) qw

      opt_flag = .false.   ! this is the simulation
      call SET_PRINTSUB(PRINTSUB) ! This is very dirty programming :)))))
#endif

C =============================================================================
C
C     INITIALIZE DATA AND MODEL
C
C =============================================================================

C     Read the *.cmp file
      call READ_CMP(MODNAM, STUBNAM, tol_daspk, tol_discret,
     1     tol_opt, zb, z0, yb, y0, ub, u0, pb, p0, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     Set derivatives of differential variables to zero
      call DCOPY(NZ, 0.d0, 0, dz, 1)

C     If file 'const_prof' exists in current directory, skip simulation and
C     copy initial estimate from *.cmp file into x
      inquire( file = 'const_prof', exist = ex )
      if( ex ) then
         call CONST_PROF(z0, dz, zb, y0, yb, u0, ub, p0, pb, x, ierr)
         write(*,75)
 75      format(/,
     1     '     File ''const_prof'' in current directory detected.',/,
     2     '     Skip simulation, and assume constant profiles.',/)
         if( ierr.ne.0 ) then
            write(*,*) 'const_prof returned IERR = ',IERR
            stop
         endif
         goto 1000
      endif

      write(*,80)
 80   format(/,'>>>>> Starting initialization of model.'/)

C     Initialize Model (e.g. for ADOL-C: create tape)
      call GET_FILENAME(MODNAM, MODNAM, ' ', filename)
!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'_model_init_' :: MODEL_INIT
      call MODEL_INIT(NZ, NY, NU, NP, TI(1), z0, dz, y0, u0, p0, f,
     1     filename)

      write(*,85)
 85   format(/,'<<<<< Done with initialization of model.'/)

#ifdef WIN32
C     Some graphics stuff
      call SET_UNIT()
#endif

C =============================================================================
C
C     DETERMINE CONSISTENT INITIAL CONDITIONS
C
C =============================================================================

C     Determine consistent initial conditions using NLEQ1S
      cnr_nleq1s = 6            ! output to screen
      tol_nleq1s = 1.d-6        ! solve relatively tight
      call INITCOND(tol_nleq1s, z0, dz, y0, u0, p0, cnr_nleq1s, ierr)
      if( ierr.ne.0 ) then
         write(*,90) ierr
 90      format(/
     1     6x,'***********************************************'/
     2     6x,'*    Problem with during NLEQ1S computation    '/
     3     6x,'*        INITCOND returned IERR = ',i4,'        *'/
     4     6x,'***********************************************'/)
         stop
      endif


C =============================================================================
C
C     SIMULATE PROCESS
C
C =============================================================================

C     Perform simulation of process using DASPK
C     Result of simulation is stored in X, with its length NX as output
C
C     Set cnr_daspk to 0, if you don't want output to file
      cnr_daspk = 6
      call DASPK_CALL(0, tol_daspk, zb, z0, dz, yb, y0, ub, u0, pb,
     1     p0, x, dummy, nx, cnr_daspk, ierr)
      if( ierr.ne.0 ) then
         if( cnr_daspk.gt.0 ) then
            write(*,100) ierr
         else
            write(*,101) ierr
         endif
 100     format(/
     1     6x,'***********************************************'/
     2     6x,'*    Problem with during DASPK simulation!    *'/
     3     6x,'*       DASPK_CALL returned IERR = ',i4,'       *'/
     4     6x,'*       Check the output file DASPK.OUT .     *'/
     4     6x,'***********************************************'/)
 101     format(/
     1     6x,'***********************************************'/
     2     6x,'*    Problem with during DASPK simulation!    *'/
     3     6x,'*       DASPK_CALL returned IERR = ',i4,'       *'/
     4     6x,'***********************************************'/)
         stop
      endif

#ifdef WIN32
!
! ------------------ Open window and display results ---------------------
!
      call open_profile_W()

! draw curves acoording to DDASSL solutions

      CALL draw_profile(0,X)

! click mouse left or right button to go ahead
      call WAIT4CLICK('LR')

! open two child windows evenly to display the outputs of calculation

C      open(outputUnit_2,file='user',
C     &            title='Collocation Information')
C      qw.x = unit_w/2.d0
C      qw.y = 0.d0
C      qw.h = unit_h - 3
C      qw.w = unit_w/2.d0
C      qw.type = QWIN$SET
C      i4 = SETWSIZEqq(OutputUnit_2,qw)
!      call clearscreen($GCLEARSCREEN)
      i4=setactiveqq(Outputunit_2)
      open(OutputUnit_1,file='user',title='IPOPT Information')
      qw.x = 1.d0
      qw.y = 0.d0
      qw.h = unit_h/2.d0 - 3
C      qw.w = unit_w/2.d0 - 4
      qw.w = unit_w - 4
      qw.type = QWIN$SET
      i4 = SETWSIZEqq(OutputUnit_1,qw)
! tile the windows opened
      i4=focusqq(outputunit_1)

      OutputUnit_COUT = OutputUnit_1
#else
      OutputUnit_COUT = 6     ! write to standard output
#endif

C =============================================================================
C
C     FIT DISCRETIZED VARIABLES
C
C =============================================================================

      inquire( file = 'skip_discret', exist = ex )
      if( .not. ex ) then
C     Refine this solution according to the discretization
         call DISCRET_IPOPT(tol_discret, nx, x, zb, yb, ub, pb, ierr)
         if( ierr.ne.0 ) then
            write(*,120) IERR
 120        format(/,'  Error in IPOPT: IERR = ',i6,'.',/)
         endif
      endif

C =============================================================================
C
C     WRITE SOLUTION
C
C =============================================================================

C     Write *.stp file containing the starting point for optimization

 1000 call WRITE_STP(.false., MODNAM, STUBNAM, x)

#ifdef WIN32
      i4=focusqq(CurveUnit)
      CALL draw_profile(1,X)
      call WAIT4CLICK('R')
      i4 = SETEXITQQ(QWIN$EXITNOPERSIST)
#endif

      return
      end
