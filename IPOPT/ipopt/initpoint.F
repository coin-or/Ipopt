C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine INITPOINT(N, X, NORIG, XORIG, IVAR, M, NIND,
     1     NLB, ILB, BNDS_L, V_L, S_L, NUB, IUB, BNDS_U, V_U, S_U, MU,
     2     CSCALE, NFIX, IFIX, KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW,
     3     IW, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     4     EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: initpoint.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Determine initial point (primal and dual) and scaling factors
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         I    INT    number of variables (without fixed vars)
CP   X        I/O   DP     I: user supplied initial point
CP                         O: corrected to keep strictly inside bounds
CP   NORIG     I    INT    number of variables in problem statement
CP                            (incl. fixed variables)
CP   XORIG    I/O   DP     I: user supplied initial point (incl. fixed vars)
CP                         O: corrected to keep strictly inside bounds
CP   IVAR      I    INT    information about partition:
CP                            X(i) corresponds to XORIG(IVAR(i))
CP   M         I    INT    number of dependent variables
CP   NIND      I    INT    number of independent variables
CP   NLB       I    INT    number of lower bounds (without fixed vars)
CP   ILB       I    INT    BNDS_L(i) is lower bound for X(ILB(i))
CP                            (Note: does not have to be ordered)
CP   BNDS_L    I    DP     lower bounds
CP   V_L       O    DP     initialized dual variables for lower bounds
CP                            (ordered as BNDS_L)
CP   S_L       O    DP     = X - BNDS_L
CP                            (ordered as BNDS_L)
CP   NUB       I    INT    number of upper bounds (without fixed vars)
CP   IUB       I    INT    BNDS_U(i) is upper bound for X(ILB(i))
CP                            (Note: does not have to be ordered)
CP   BNDS_U    I    DP     upper bounds
CP   V_U       O    DP     initialized dual variables for lower bounds
CP                            (ordered as BNDS_U)
CP   S_U       O    DP     = X - BNDS_U
CP                            (ordered as BNDS_U)
CP   MU        O    DP     initial value of barrier parameter
CP   CSCALE    O    DP     scaling factors (see get_scale)
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS        S    DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS        S    INT    INT storage space (all!)
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective fuction
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      include 'IPOPT.INC'

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      double precision X(N)
      integer NORIG
      double precision XORIG(NORIG)
      integer IVAR(N)
      integer M
      integer NIND
      integer NLB
      integer ILB(NLB)
      double precision BNDS_L(NLB)
      double precision V_L(NLB)
      double precision S_L(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_U(NUB)
      double precision V_U(NUB)
      double precision S_U(NUB)
      double precision MU
      double precision CSCALE(*)
      integer NFIX
      integer IFIX(NFIX)
      integer KCONSTR(6)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision DNRM2, DASUM, D1MACH
      integer IDAMAX

      integer p_rwend, p_iwend, p_ilbinv, p_iubinv, p_g, p_c
      integer i, ibl, ibu, ii, idummy
      double precision bl, bu, diff, xtmp, vinit, macheps, fscale_orig
      double precision dummy, xi, xscale
      logical upp, low
      character*80 line(3)

C
C     QINIT   =  1: use MU=QMUINIT, v = 1d0
C               -1: use MU=QMUINIT, v = MU/s
C                0: use MU=QMUINIT, take v's as given by user
C                2: v = max(1d0, ||G(x0)||_2), MU = max(1, (S^T V)/(NLB+NUB))
C
C     QSCALE  =  0: use QFSCALE as in PARAMS.DAT, don't scale constraints
C                1: use Yamashita heuristic to overwrite QFSCALE and to obtain
C                   scalar scaling factor CSCALE for ALL constraints
C                2: read vector of scaling factor for constraints from file
C                   CSCALE.DAT
C
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
      macheps = D1MACH(4)
      p_rwend = 0
      p_iwend = 0
C
C     Output
C
      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
         write(line,700)
 700     format(/,'Information about initial point:',/)
         call C_OUT(1,2,3,line)
      endif
C
C     Compute initial values of slacks; if original point too close to or
C        violating bounds, X-S not 0.
C

C     Determine inverse of ILB and IUB
      p_ilbinv = p_iwend
      p_iubinv = p_ilbinv + N
      p_iwend  = p_iubinv + N
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      do i = 1, N
         IW(p_ilbinv+i) = 0
         IW(p_iubinv+i) = 0
      enddo
      do i = 1, NLB
         IW(p_ilbinv+ILB(i)) = i
      enddo
      do i = 1, NUB
         IW(p_iubinv+IUB(i)) = i
      enddo

C     Loop over all variables to determine if they have lower and/or upper
C        bounds.  If so, initialize slacks

      do i = 1, N
         ibl = IW(p_ilbinv+i)
         if( ibl.ne.0 ) then
            low = .true.
CTODO Decide if there should be an option to move bounds slightly
            if( QMOVEBOUNDS.gt.0.d0 ) then
               BNDS_L(ibl) = BNDS_L(ibl) -
     1              QMOVEBOUNDS*max(1.d0,dabs(BNDS_L(ibl)))
            endif
            bl = BNDS_L(ibl)
         else
            low = .false.
         endif
         ibu = IW(p_iubinv+i)
         if( ibu.ne.0 ) then
            upp = .true.
            if( QMOVEBOUNDS.gt.0.d0 ) then
               BNDS_U(ibu) = BNDS_U(ibu) +
     1              QMOVEBOUNDS*max(1.d0,dabs(BNDS_U(ibu)))
            endif
            bu = BNDS_U(ibu)
         else
            upp = .false.
         endif

         if( QSCALE.ge.3 ) then
            xscale = CSCALE(M+i)
         else
            xscale = 1.d0
         endif

         if( low ) then
            if( upp ) then
               diff = bu - bl
               if( diff.le.0.d0 ) then
                  write(line(1),*) 'Bounds ',i,' are inconsistent:'
                  write(line(2),*) 'BL = ',bl,', BU = ',bu
                  call C_OUT(2,0,2,line)
                  IERR = 7
                  goto 9999
               endif
               xi = X(i)
C#define OLD
#ifdef OLD
               xtmp = dmax1(dmin1(bl + QFRAC*diff,bl + QTHETA),
     1                      dmin1(dmax1(bu - QFRAC*diff, bu - QTHETA),
     1                            X(i)))
#else
               xtmp = dmax1(dmin1(bl + QFRAC*diff,
     1              bl + QTHETA*max(1.d0,abs(bl))),
     1              dmin1(dmax1(bu - QFRAC*diff,
     1              bu - QTHETA*max(1.d0,abs(bu))), X(i)))
#endif
               X(i) = xtmp
               if( xi.ne.xtmp ) then
                  ii = IVAR(i)
                  if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
                     write(line,710) ii,XORIG(ii),xtmp*xscale
 710                 format('  XORIG(',i7,') corrected from ',
     1                      d12.5,' to ',d12.5)
                     call C_OUT(1,2,1,line)
                  endif
                  XORIG(ii) = X(i)*xscale
               endif
               S_L(ibl) = xtmp - bl
               S_U(ibu) = bu - xtmp
            else
               xi = X(i)
#ifdef OLD
               xtmp = dmax1(bl + QTHETA, X(i))
#else
               xtmp = dmax1(bl + QTHETA*max(1.d0,abs(bl)), X(i))
#endif
               X(i) = xtmp
               if( xi.ne.xtmp ) then
                  ii = IVAR(i)
                  if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
                     write(line,710) ii,XORIG(ii),xtmp*xscale
                     call C_OUT(1,2,1,line)
                  endif
                  XORIG(ii) = X(i)*xscale
               endif
               S_L(ibl) = xtmp - bl
            endif
         else
            if( upp ) then
               xi = X(i)
#ifdef OLD
               xtmp = dmin1(bu - QTHETA, X(i))
#else
               xtmp = dmin1(bu - QTHETA*max(1.d0,abs(bu)), X(i))
#endif
               X(i) = xtmp
               if( xi.ne.xtmp ) then
                  ii = IVAR(i)
                  if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
                     write(line,710) ii,XORIG(ii),xtmp*xscale
                     call C_OUT(1,2,1,line)
                  endif
                  XORIG(ii) = X(i)*xscale
               endif
               S_U(ibu) = bu - xtmp
            endif
         endif
      enddo
 190  continue
      p_iwend = p_ilbinv

      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then
         call C_OUT(1,2,1,'  All slacks satisfy slack equation.')
      endif
C
C     Initialize barrier parameter and dual variables
C
      if( abs(QINIT).le.1 ) then
C
C     Initialize barrier parameter
C
         MU = QMUINIT
C
C     Initialize slack multipliers (now as MU*S^{-1}e )
C
         if( QINIT.eq.-1 ) then
            do i = 1, NLB
               V_L(i) = MU/S_L(i)
            enddo
            do i = 1, NUB
               V_U(i) = MU/S_U(i)
            enddo
         elseif( QINIT.eq.1 ) then
            vinit = 1.d0
            call DCOPY(NLB, vinit, 0, V_L, 1)
            call DCOPY(NUB, vinit, 0, V_U, 1)
         else                   ! warm start
            do i = 1, NLB
               V_L(i) = dmax1( V_L(i), QTHETA )
            enddo
            do i = 1, NUB
               V_U(i) = dmax1( V_U(i), QTHETA )
            enddo
         endif
C
      elseif( QINIT.eq.2 ) then
C
C     Initialize dual variables based on norm of objective gradient
C
         if( NLB+NUB.eq.0 ) then
            MU = QMUINIT
         else
            p_g     = p_rwend
CTODO            p_rwend = p_g + N
            p_rwend = p_g + N + M
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, 0, idummy,
     1           dummy, 0, idummy, dummy, 0.d0, RW(p_g+1),
     1           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     1           IW(p_iwend+1), IERR, EV_G, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*) 'get_scale: get_g returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C     IF LEAVE HERE: CHECK NAN INF AND REMOVE FROM MAINLOOP?
      call GET_C(0, N, NIND, X, IVAR, NORIG, XORIG, M,
     1     CSCALE, RW(p_g+N+1), KCONSTR, LRS, RS, LIS, IS,
     2     LRW-p_rwend, RW(p_rwend+1),
     4     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         write(line,*)
     1        'mainloop: Error: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      elseif( IERR.ne.0 ) then
         write(line,*)
     1        'mainloop: Warning: get_c returns IERR = ',IERR
         call C_OUT(2,0,1,line)
      endif
CTODO decide which norm (l1 or l2) to use...
            vinit = DNRM2(N, RW(p_g+1), 1)/sqrt(dble(NLB+NUB))
            i = IDAMAX(N, RW(p_g+1), 1)
            vinit = abs(RW(p_g+i))
            vinit = DASUM(N+M, RW(p_g+1), 1)/dble(NLB+NUB+M)
c            vinit = DASUM(N, RW(p_g+1), 1)
            p_rwend = p_g
C         vinit = dmax1( vinit, 1d0 )
            vinit = min(1.d3, max( vinit, 1d-3 ))
            write(line,*) 'initpoint: vinit = ',vinit
            call C_OUT(2,2,1,line)
            call DCOPY(NLB, vinit, 0, V_L, 1)
            call DCOPY(NUB, vinit, 0, V_U, 1)
C
C     Compute initial value of barrier parameter
C
            MU = 0.d0
            do i = 1, NLB
               MU = MU + V_L(i)*S_L(i)
            enddo
            do i = 1, NUB
               MU = MU + V_U(i)*S_U(i)
            enddo
C            MU = dmax1( QMUINIT, MU/(NLB+NUB) )
            MU = min(1.d3, max( 1.d-1*MU/dble(NLB+NUB), 1.d-3 ))
c            MU = QMUINIT
C         do i = 1, NLB
C            V_L(i) = min( 1.d3, max( 1.d-3, MU/S_L(i) ))
C         enddo
C         do i = 1, NUB
C            V_U(i) = min( 1.d3, max( 1.d-3, MU/S_U(i) ))
C         enddo
            write(line,*) 'initpoint: MUinit  = ',MU
            call C_OUT(2,2,1,line)
         endif
      elseif( QINIT.eq.3 ) then
         MU = QMUINIT
         do i = 1, NLB
            V_L(i) = min( 1.d3, max( 1.d-3, MU/S_L(i) ))
         enddo
         do i = 1, NUB
            V_U(i) = min( 1.d3, max( 1.d-3, MU/S_U(i) ))
         enddo
      elseif( QINIT.eq.4 .or. QINIT.eq.5 ) then
         if( NLB+NUB.eq.0 ) then
            MU = QMUINIT
         else
            vinit = 1.d0
            call DCOPY(NLB, vinit, 0, V_L, 1)
            call DCOPY(NUB, vinit, 0, V_U, 1)
            MU = 0.d0
            do i = 1, NLB
               MU = MU + V_L(i)*S_L(i)
            enddo
            do i = 1, NUB
               MU = MU + V_U(i)*S_U(i)
            enddo
C            MU = min(1.d6, max( MU/dble(NLB+NUB), 1.d-4 ))
CTODO Don't need lower bound if vinit is bounded, since minimal distance to
C     boundary for slacks imply lower bound
            MU = min(1.d6, MU/dble(NLB+NUB))
            write(line,*) 'initpoint: MUinit  = ',MU
            call C_OUT(2,2,1,line)
            if( QINIT.eq.5 ) then
               do i = 1, NLB
                  V_L(i) = min( 1.d3, max( 1.d-3, MU/S_L(i) ))
               enddo
               do i = 1, NUB
                  V_U(i) = min( 1.d3, max( 1.d-3, MU/S_U(i) ))
               enddo
            endif
         endif
      else
         write(line,*) 'initpoint: QINIT = ',QINIT
         call C_OUT(2,0,1,line)
         IERR = 4
         goto 9999
      endif
C
C     Make sure that dual variables don't differ too much from MU*X{-1}*e
C
      call CORRECT_V(NLB, S_L, V_L, NUB, S_U, V_U, MU)
CC
CC     If we are dealing with a square problem, we choose MU to be less than
CC     error tolerance, so that there are no initial dummy iterations
CC
C      if( NIND.eq.0 .or. (NLB.eq.0 .and. NUB.eq.0) .and. QCG.eq.0 ) then
C         MU = 0.5d0*QTOL/QMUERRFAC
C      endif

C     Set common block variable for AMPL interface
      AMPLMU = MU

 9999 continue
      return
      end

      subroutine CORRECT_V(NLB, S_L, V_L, NUB, S_U, V_U, MU)
C
C     Subroutine correcting V is necessary to ensure that the primal-dual
C     Hessian Sigma does not deviate too much from the primal-only Hessian
C
      implicit none
      include 'IPOPT.INC'
      integer NLB, NUB
      double precision S_L(NLB), V_L(NLB), S_U(NUB), V_U(NUB), MU

      double precision tmp
      character*80 line
      integer i
C
C     Make sure that dual variables don't differ too much from MU*X{-1}*e
C
      do i = 1, NLB
         tmp = V_L(i)
         V_L(i) = dmax1(V_L(i), MU/QVCORRECTFACT/S_L(i))
         V_L(i) = dmin1(V_L(i), MU*QVCORRECTFACT/S_L(i))
         if( tmp.ne.V_L(i).and.QPRINT.gt.2 ) then
            write(line,*) 'correct_v: V_L(',i,') corrected from'
            call C_OUT(1,0,1,line)
            write(line,*) '           ',tmp,' to ',V_L(i)
            call C_OUT(1,0,1,line)
         endif
      enddo
      do i = 1, NUB
         tmp = V_U(i)
         V_U(i) = dmax1(V_U(i), MU/QVCORRECTFACT/S_U(i))
         V_U(i) = dmin1(V_U(i), MU*QVCORRECTFACT/S_U(i))
         if( tmp.ne.V_U(i).and.QPRINT.gt.2 ) then
            write(line,*) 'correct_v: V_U(',i,') corrected from'
            call C_OUT(1,0,1,line)
            write(line,*) '           ',tmp,' to ',V_U(i)
            call C_OUT(1,0,1,line)
         endif
      enddo
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine INITPOINT_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1

      LRW = 0
      LIW = 2*N

      if( QINIT.eq.2 ) then
         call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
         LRW = max(LRW, N+lrw1)
         LIW = max(LIW, liw1)
      endif

      return
      end
