C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine ITER_OUT(ITER, MU, ERR, C_ERR, ERR_RDGRD, ERR_CNSTR,
     1     ERR_CMPL, F, M, C, CNRM0, YPY, NIND, PZ, C_SKIP,
     1     ALFAV, ALFATAU, ALFAX, C_ALPHA, C_WATCH, LS_COUNT,
     3     NU, IEIGS, CONDC, REGU,
     2     RV, RG, RGB, B, SIGMA_L, SIGMA_U, W, ZPZ, DX,
     3     X, S_L, S_U, V_L, V_U, N, NLB, ILB, NUB, IUB, IVAR,
     4     LAM, XORIG)
C
C*******************************************************************************
C
C    $Id: iter_out.F 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Output of information of iteration
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
CP   ITER      I    INT    iteration counter
CP   MU        I    DP     barrier parameter
CP   ERR       I    DP     KKT error (if QERROR < 0 scaled)
CP   C_ERR     I    C*1    info about error
CP   ERR_RDGRD I    DP     unscaled/scaled max norm of dual infeasibility
CP                            (red. gradient in case of reduced space approach)
CP   ERR_CNSTR I    DP     unscaled/scaled max norm of infeasibility
CP   ERR_CMPL  I    DP     unscaled/scaled max norm for relaxed complementarity
CP   F         I    DP     value of objective function (at new point after step)
CP   M         I    INT    number of constraints
CP   C         I    DP     values of constraints (at new point after step)
CP   CNRM0     I    DP     2-norm of constraints at current point
CP   YPY       I    DP     range space step
CP   NIND      I    INT    number of inpendent variables
CP   PZ        I    DP     null space step
CP   C_SKIP    I    C*1    info about Quasi-Newton update
CP   ALFAV     I    DP     step size for dual variables
CP   ALFATAU   I    DP     step size for primal variables before line search
CP   ALFAX     I    DP     step size for primal variables after line search
CP   C_ALPHA   I    C*1    info about first alpha
CP   C_WATCH   I    C*1    info about watchdog
CP   LS_COUNT  I    INT    number of trial steps in line search
CP   NU        I    DP     value of penalty parameter (for QQUASI = -2,-1,1,2,3)
CP   IEIGS     I    INT    info about corrections if red. Hessian indefinite /
CP                            cg iterations
CP   CONDC     I    DP     estimate of condition number
CP   REGU      I    DP     regularization perturbation of Hessian
CP   RV        I    DP     reduced multipliers
CP   RG        I    DP     reduced gradient of objective function
CP   RGB       I    DP     reduced primal-dual gradient of barrier function
CP   B         I    DP     reduced Hessian (estimate) of original Lagrangian
CP   SIGMA_L   I    DP     primal-dual Hessian of barrier term for lower bounds
CP   SIGMA_U   I    DP     primal-dual Hessian of barrier term for upper bounds
CP   W         I    DP     overall reduced Hessian
CP   ZPZ       I    DP     overall null space step
CP   DX        I    DP     overall primal step
CP   X         I    DP     iterate at new point
CP   S_L       I    DP     slacks for lower bounds at new point
CP   S_U       I    DP     slacks for upper bounds at new point
CP   V_L       I    DP     multipliers for lower bounds at new point
CP   V_U       I    DP     multipliers for upper bounds at new point
CP   N         I    INT    number of (free) variables
CP   NLB       I    INT    number of lower bounds
CP   ILB       I    INT    indices for lower bounds
CP   NUB       I    INT    number of upper bounds
CP   IUB       I    INT    indices for upper bounds
CP   IVAR      I    INT    parition information
CP   LAM       I    DP     value of multipliers for equality constraints
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
CS    DNRM2
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
      integer ITER
      double precision MU
      double precision ERR
      character*1 C_ERR
      double precision ERR_RDGRD
      double precision ERR_CNSTR
      double precision ERR_CMPL
      double precision F
      integer M
      double precision C(M)
      double precision CNRM0
      integer N
      double precision YPY(N)
      integer NIND
      double precision PZ(NIND)
      character*1 C_SKIP
      double precision ALFAV
      double precision ALFATAU
      double precision ALFAX
      character*1 C_ALPHA
      character*1 C_WATCH
      integer LS_COUNT
      double precision NU
      integer IEIGS
      double precision CONDC
      double precision REGU
      double precision RV(NIND)
      double precision RG(NIND)
      double precision RGB(NIND)
      double precision B((NIND*(NIND+1))/2)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
      double precision W((NIND*(NIND+1))/2)
      double precision ZPZ(M)
      double precision DX(N)
      double precision X(N)
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      integer IVAR(N)
      double precision LAM(M)
      double precision XORIG(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision DNRM2, DDOT
      double precision pynrm, zpznrm, dnrm, cputim
      integer i
      character*150 line(4)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( QFULL.eq.0 .and. ITER.ne.0 ) then
         pynrm = DNRM2(N, YPY, 1)
         zpznrm = DDOT(NIND, PZ, 1, PZ, 1)
C      if( QFULL.ne.1 ) then
C         zpznrm = zpznrm + DDOT(M, ZPZ, 1, ZPZ, 1)
C      endif
         zpznrm = DSQRT(zpznrm)
      else
         pynrm = 0.d0
         zpznrm = 0.d0
      endif
      if( ITER.ne.0 ) then
         dnrm = DNRM2(N, DX, 1)
      else
         dnrm = 0.d0
      endif

      if( QFULL.ne.0 .or. QCG.eq.0 ) then
         if(QOUTPUT.eq.0) then
            write(line,201)
         else
            write(line,101)
         endif
      else
         if(QOUTPUT.eq.0) then
            write(line,202)
         else
            write(line,102)
         endif
      endif
      if( mod(ITER,10).eq.0 ) then
         if( ITER.eq.0 .and. LASTPRINT.eq.0 ) then
            call C_OUT(1,2,2,line)
         else
            call C_OUT(2,1,2,line)
            LASTPRINT = ITER
         endif
      else
         call C_OUT(1,2,2,line)
      endif

      call TIMER(cputim)
      if( QOUTPUT.eq.0 ) then
         write(line,211)
     1        ITER, ERR, C_ERR, MU, CNRM0, dnrm, ALFAX, C_WATCH,
     1        LS_COUNT, F/QFSCALE, REGU
         call C_OUT(0,1,1,line)
         write(line,111)
     1        ITER, ERR, C_ERR, MU, CNRM0, pynrm, zpznrm,
     1        C_SKIP, dnrm, ALFAV, C_ALPHA, ALFAX, C_WATCH, NU,
     1        LS_COUNT, F/QFSCALE, IEIGS, REGU, cputim
         call C_OUT(1,1,1,line)
      else
         write(line,111)
     1        ITER, ERR, C_ERR, MU, CNRM0, pynrm, zpznrm,
     1        C_SKIP, dnrm, ALFAV, C_ALPHA, ALFAX, C_WATCH, NU,
     1        LS_COUNT, F/QFSCALE, IEIGS, REGU, cputim
         call C_OUT(2,1,1,line)
      endif

 201  format(/,'ITER     ERR       MU      ||C||    ||D||   ',
     1     'ALFA(X) #LS        F         Regu')

 101  format(/,'ITER     ERR       MU      ||C||   ||',
     1           'YPY||  ||PZ||     ||D||   ALFA(V)',
     2           '   ALFA(X)  ',
     2           '   NU   #LS       F      #cor   Regu    CPU(s)')

 202  format(/,'ITER     ERR        MU       ||C||     ||D||     ',
     1     'ALFA(X)  #LS        F      #cg')

 102  format(/,'ITER     ERR       MU      ||C||   ||',
     1           'YPY||  ||PZ||     ||D||   ALFA(V)',
     2           '   ALFA(X)  ',
     2           '   NU   #LS       F       #cg   Regu    CPU(s)')

 211  format(i5,1x,d8.3,a1,1x,d8.3,1x,d8.3,1x,d8.3,1x,d8.3,a1,
     1       i3,1x,d14.8,1x,d8.3)

 111  format(i5,1x,d8.3,a1,1x,d8.3,1x,d8.3,1x,d8.3,1x,d8.3,a1,1x,d8.3,
     1       1x,d8.3,a1,1x,d8.3,a1,1x,d8.3,i3,1x,d14.8,i3,1x,d8.3,1x,
     2       d8.3)

C
C     If wished write even more info into a file
C
      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then

         write(line,700) ITER
 700     format(/,80('-'),/,20(' '),
     1          'Detailed Information for Iteration ',
     1          i7,/,80('-'))
         call C_OUT(1,0,4,line)

         write(line,705) ITER,F/QFSCALE
 705     format(/,'  Objective function (Iter ',i7,'):',//,
     1        ' F = ',d25.17)
         call C_OUT(1,0,4,line)

         if( QPRINT.ge.9 .and. QFULL.eq.0 ) then
            write(line,710) ITER
 710        format(/,'  Reduced bound multipliers (Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NIND
               write(line,711) i,RV(i)
 711           format(' RV(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.8 .and. (QFULL.eq.0 .or. QQUASI.ne.0) ) then
            write(line,720) ITER
 720        format(/,'  Reduced gradient of objective function ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NIND
               write(line,721) i,RG(i)
 721           format(' RG(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.8 .and. (QFULL.eq.0 .or. QQUASI.ne.0) )then
            write(line,730) ITER
 730        format(/,'  Reduced gradient of barrier function ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NIND
               write(line,731) i,RGB(i)
 731           format(' RGB(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.7 .and. ITER.gt.0 .and. QQUASI.ne.0 .and.
     1        abs(QQUASI).le.0)then
            write(line,740) ITER
 740        format(/,'  Quasi Newton Matrix ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, (NIND*(NIND+1))/2
               write(line,741) i,B(i)
 741           format(' B(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.8 )then
            write(line,750) ITER
 750        format(/,'  Sigmas for lower and upper bounds ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NLB
               write(line,751) IVAR(ILB(i)),SIGMA_L(i)
 751           format(' SIGMA_L(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
            call C_OUT(1,0,0,line)
            do i = 1, NUB
               write(line,752) IVAR(IUB(i)),SIGMA_U(i)
 752           format(' SIGMA_U(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.7 .and. ITER.gt.0 .and. QFULL.eq.0 .and.
     1        QQUASI.le.5 )then
            write(line,760) ITER
 760        format(/,'  reduced Hessian ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, (NIND*(NIND+1))/2
               write(line,761) i,W(i)
 761           format(' W(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 .and. QFULL.eq.0 ) then
            write(line,770) ITER
 770        format(/,'  range space step ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, N
               write(line,771) i,YPY(i)
 771           format(' YPY(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 .and. QFULL.eq.0 )then
            write(line,780) ITER
 780        format(/,'  null space step (indepentend vars) ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NIND
               write(line,781) i,PZ(i)
 781           format(' PZ(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 .and. QFULL.eq.0 ) then
            write(line,790) ITER
 790        format(/,'  null space step (depentend vars) ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, M
               write(line,791) i,ZPZ(i)
 791           format(' ZPZ(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 )then
            write(line,795) ITER
 795        format(/,'  overall step ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, N
               write(line,796) IVAR(i),DX(i)
 796           format(' DX(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,800) ITER
 800        format(/,'  new primal variables ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, N
               write(line,801) IVAR(i),X(i)
 801           format(' X(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,810) ITER
 810        format(/,'  new slack variables ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NLB
               write(line,811) IVAR(ILB(i)),S_L(i)
 811           format(' S_L(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
            call C_OUT(1,0,0,line)
            do i = 1, NUB
               write(line,812) IVAR(IUB(i)),S_U(i)
 812           format(' S_U(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,820) ITER
 820        format(/,'  new dual variables ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NLB
               write(line,821) IVAR(ILB(i)),V_L(i)
 821           format(' V_L(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
            call C_OUT(1,0,0,line)
            do i = 1, NUB
               write(line,822) IVAR(IUB(i)),V_U(i)
 822           format(' V_U(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.5 .and. QLAMBDA.ne.0 )then
            write(line,825) ITER
 825        format(/,'  new multipliers for equalities ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, M
               write(line,826) i,LAM(i)
 826           format(' LAM(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,827) ITER
 827        format(/,'  new constraint values ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, M
               write(line,828) i,C(i)
 828           format(' C(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif


         if( QPRINT.ge.7 )then
            write(line,830) ITER
 830        format(/,'  complementarity: ',
     1             '(Iter ',i7,'):',/)
            call C_OUT(1,0,3,line)
            do i = 1, NLB
               write(line,831) IVAR(ILB(i)),IVAR(ILB(i)),V_L(i)*S_L(i)
 831           format(' V_L(',i7,')*S_L(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
            call C_OUT(1,0,0,line)
            do i = 1, NUB
               write(line,832) IVAR(IUB(i)),IVAR(IUB(i)),V_U(i)*S_U(i)
 832           format(' V_U(',i7,')*S_U(',i7,') = ',d25.17)
               call C_OUT(1,0,1,line)
            enddo
         endif

      endif

      return
      end
