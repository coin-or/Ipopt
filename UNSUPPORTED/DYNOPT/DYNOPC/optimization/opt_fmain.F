C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C =============================================================================
C $Id: opt_fmain.F 531 2004-03-11 01:31:07Z andreasw $
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
      external DRAW_ITER_INFO
C !DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
#endif

C     Bounds and initial values for model variables read from .cmp file
      double precision zb(2,NZMAX), z0(NZMAX)
      double precision yb(2,NYMAX), y0(NYMAX)
      double precision ub(2,NUMAX), u0(NUMAX)
      double precision pb(2,NPMAX), p0(NPMAX)

C     Derivatives for differential variables
      double precision dz(NZMAX)

C     Dummy residual for initialization of MODEL (MODEL_INIT)
      double precision f(NZMAX+NYMAX)

C     Tolerances read from .cmp file
      double precision tol_daspk, tol_discret, tol_opt

C     Vector X for discretization
      double precision x(NXMAX)

C     Time discretization from starting point file
      double precision ti_stp(NELEMAX+1)
      integer nele_stp, ncol_stp

C     Defitinions for IPOPT
      integer n, m
      integer nlb, nub
      integer ilb(NXMAX), iub(NXMAX)
      double precision bnds_l(NXMAX), bnds_u(NXMAX)
      double precision v_l(NXMAX), v_u(NXMAX), lam(NXMAX)
      integer iter
      integer narg
      double precision args(50)
      character*20 cargs(50)

      double precision rw(LRW_IPOPT)
      integer          iw(LIW_IPOPT)
C      integer PRW, PIW, MYALLOC

      integer ierr, nlines, ndegu, i
      character*80 lines(7)
      character*256 fname
      logical warm, ex, new
      character*1 answer, s
      double precision d, times, timef
      integer NoCorba

      integer READ_ENTRY, CALC_NX

#ifdef WIN32
      integer(4) i4
      type(qwinfo) qw

      call SET_PRINTSUB(PRINTSUB)
      call SET_DRAWSUB(DRAW_ITER_INFO)
#endif

C =============================================================================
C
C     Read Data Files
C
C =============================================================================

      write(*,100)
 100  format(/
     1     '>>>>> Starting reading input data.')

C     Read the *.cmp file
      call READ_CMP(MODNAM, STUBNAM, tol_daspk, tol_discret,
     1     tol_opt, zb, z0, yb, y0, ub, u0, pb, p0, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     Check, whether there is a *.sol file.  If so, ask if user wants
C     to use this as starting point ('warm start')
      call GET_FILENAME(MODNAM, STUBNAM, 'sol', fname)
      inquire(file=fname, exist=warm)
      if( warm ) then
#ifdef WIN32
         write(*,102)
 102     format(/,'  Do you want a warm start? (y/n): '\)
#else
         call CPRINTF('\n%s','  Do you want a warm start? (y/n): ')
#endif
         read(*,*) answer
C         answer = 'n'
         warm = (answer.eq.'Y' .or. answer.eq.'y')
      endif

C     Read the starting point
      call READ_STP(warm, MODNAM, STUBNAM, x, ti_stp,
     1     nele_stp, ncol_stp, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         write(*,*)
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     If necessary, fit starting point to a new TI grid
      call UPDATE_STP(x, ti_stp, nele_stp, ncol_stp, ierr)

C     Read information on the objective function
      call READ_OJF(MODNAM, STUBNAM, ndegu, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         write(*,*)
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     Read information on additional bound constraints
      call READ_ADB(MODNAM, STUBNAM, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         write(*,*)
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C     Determine element, on which factorization is to be based
      call GET_FILENAME(MODNAM, STUBNAM, 'badele', fname)
      inquire(file=fname, exist=ex)
      if( ex ) then
         open(1,file=fname,status='old')
         new = .true.
         if( READ_ENTRY(1, new, 'i', d, CRIT_ELE, s).ne.0 ) then
            write(*,105)
            stop
         endif
         if( CRIT_ELE.le.0 .or. CRIT_ELE.gt.NELE ) then
            write(*,105)
            stop
         endif
 105     format(/,6x,'Error reading "*.badele" file.  Abort.')
         close(1,status='keep')
      else
         CRIT_ELE = NELE/2 + 1
      endif

      write(*,110)
 110  format(/
     1     '<<<<< Done reading input data.'/)

C =============================================================================
C
C     Initialize Data for IPOPT
C
C =============================================================================

C     Set initial conditions
      call DCOPY(NZ, z0, 1, ZINIT, 1)

C     Possibly adjust bounds in order to ensure (well heuristic) a strict
C     interior
      call PUSH_BND(tol_opt, zb, yb, x)

C     Set bounds (NLB, NUB, ILB, IUB, BNDS_L, BNDS_U)
      call SET_BNDS(zb, yb, ub, pb, nlb, nub, ilb, iub,
     1     bnds_l, bnds_u)

C
CDELETEME
      if( .false. ) then
      write(*,*) 'BOUNDS FOR TE included!'
      do i = 1, NZ
         nlb = nlb + 1
         bnds_l(nlb) = - 1e-6 *max(1,abs(ZINIT(i)))
         ilb(nlb) = NZ + (NELE-1)*(NCOL*(NZ+NY)+NZ)
     1        + (NCOL-1)*(NZ+NY) + i
         nub = nub + 1
         bnds_u(nub) = + 1e-6 *max(1,abs(ZINIT(i)))
         iub(nub) = ilb(nlb)
      enddo
      endif

C     Set additional constraints to enforce controls of degree ndegu
      call ADDCON_INIT(NZ, NY, NU, NCOL, ndegu, NAC)

C     Non of the controls and profiles are fixed
      NU_PROF = 0
      NP_FIX  = 0
      do i = 1, NU
         IU_PROF(i) = 0
      enddo
      do i = 1, NP
         IP_FIX(i) = 0
      enddo

C =============================================================================
C
C     Initialize Model
C
C =============================================================================

      write(*,150)
 150  format(/
     1     '>>>>> Starting initializing model.'/)

C     Initialize Model (e.g. for ADOL-C: create tape)
      call GET_FILENAME(MODNAM, MODNAM, ' ', fname)
C !DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'_model_init_' :: MODEL_INIT
      call MODEL_INIT(NZ, NY, NU, NP, TI(1), z0, dz, y0, u0, p0, f,
     1     fname)

      write(*,160)
 160  format(/
     1     '<<<<< Done initializing model.')

#ifdef WIN32
C =============================================================================
C
C     Initialize Graphic
C
C =============================================================================

! open a frame window and set it to maximum
      opt_flag = .true.

      qw.type=QWIN$MAX
      I4=setwsizeqq(QWIN$FRAMEWINDOW,qw)
      CALL SET_UNIT_OPT()
!      open child windows and calculate unit of width and height based on character
      CALL open_profile_W()
      call open_alfa_w()
      OPT_EPS = tol_opt         !  in /GRAPH/
      call open_rginf_w() 

      i4 = setactiveqq(CurveUnit)
      call DRAW_PROFILE(1, x)
#endif

C =============================================================================
C
C     Call IPOPT
C
C =============================================================================

C     Set IPOPT parameters
      cargs(1) = 'ifull'
      args(1) = 0.d0
      cargs(2) = 'ihessvect'
      args(2) = 2.d0
      cargs(3) = 'dtol'
      args(3) = tol_opt
      narg = 3
      call INITPARAMS(ierr, narg, args, cargs)
      if( ierr.ne.0 ) then
         write(*,200) IERR
 200     format(/,6x,'Error initializing IPOPT parameters. IERR =',i6)
         return
      endif

      write(*,210)
 210  format(/
     1     '>>>>> Calling IPOPT to perform the optimization.'/)

      n = CALC_NX(NZ, NY, NU, NP, NELE, NCOL)
      m = n - NELE*(NCOL*NU - NAC) - NP

C     Set counters for model calls to zero
      NO_MODEL_F_CALL = 0
      NO_MODEL_DF_CALL = 0

      call TIMER(times)

      call IPOPT(n, x, m, nlb, ilb, bnds_l, nub, iub, bnds_u, v_l, v_u,
     1     lam, LRW_IPOPT, rw, LIW_IPOPT, iw, iter, IERR)

C      PRW = MYALLOC(8*LRW_IPOPT)
C      PIW = MYALLOC(4*LIW_IPOPT)
C      call IPOPT(n, x, m, nlb, ilb, bnds_l, nub, iub, bnds_u, v_l, v_u,
C     1     lam, LRW_IPOPT, %VAL(PRW), LIW_IPOPT, %VAL(PIW), iter, IERR)

      call TIMER(timef)

      write(*,215) timef-times

 215  format(/,6x,'CPU seconds spend in IPOPT = ',d11.4)

      write(*,250) NO_MODEL_F_CALL, NO_MODEL_DF_CALL
 250  format(/,6x,'Number of model function evaluations = ',i6,
     1       /,6x,'Number of model Jacobian evaluations = ',i6)
c      Call Get_CorbaCalls(NoCorba)
c      write(*,255) NoCorba
c 255  format(/,6x,'               Number of CORBA Calls = ',i6)

      if( IERR.eq.0 ) then
         write(*,220)
 220     format(/
     1        '<<<<< IPOPT successfully optimized the process.'/)
      else
         write(*,230) IERR
 230     format(/
     1        6x,'IPOPT did not converge  --- IERR = ',i6,'.'/)
      endif

#ifdef WIN32
      call open_report_w()
      call draw_profile(1,x)
#endif

DTODO:

C =============================================================================
C
C     Call DASPK to see if discretization is fine enough
C     (At this point, we assume controls constant in each element)
C
C =============================================================================

C      cnr_daspk = 6
C      tol_daspk = 1d-8          ! choose very fine tolerance for now
C      call DASPK_CALL(1, tol_daspk, zb, z0, dz, yb, y0, ub, u0,
C     1     pb, p0, x2, x, nx, cnr_daspk, ierr)
C      if( ierr.ne.0 ) then
C         write(*,250) IERR
C 250     format(/,6x,'Simulation with DASPK not successful',/,
C     1        'DASPK_CALL returns IERR = ',i6,/)
C      endif


      write(*,300)
 300  format(/
     1     '>>>>> Starting writing output data.')

C     Store solution of optimization
      call WRITE_STP(.true., MODNAM, STUBNAM, x )

      write(*,310)
 310  format(/
     1     '<<<<< Done writing output data.'/)

#ifdef WIN32
      call wait4click('LR')
#endif

      return
      end
