C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C*******************************************************************************
C
      subroutine ITER_OUT(ITER, MU, ERR, C_ERR, ERR_RDGRD, ERR_CNSTR,
     1     ERR_CMPL, F, M, C, CNRM0, YPY, NIND, PZ, C_SKIP,
     1     ALFAV, ALFATAU, ALFAX, C_ALPHA, C_WATCH, LS_COUNT,
     3     NUH, IEIGS, CONDC, LSLACKS,
     2     RV, RG, RGB, B, SIGMA_L, SIGMA_U, W, ZPZ, DX,
     3     X, S_L, S_U, V_L, V_U, N, NLB, ILB, NUB, IUB, IVAR,
     4     LAM, XORIG)
C
C*******************************************************************************
C
C    $Id: iter_out_dynopc.F 531 2004-03-11 01:31:07Z andreasw $
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
CA    Andreas Waechter      02/22/99
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
CP   NUH       I    DP     value of penalty parameter (for QQUASI = -2,-1,1,2,3)
CP   IEIGS     I    INT    info about corrections if red. Hessian indefinite /
CP                            cg iterations
CP   CONDC     I    DP     estimate of condition number
CP   LSLACKS   I    LOG    =.true.: slacks violate linear slack constraint
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
#ifdef WIN32
      include 'DYNAUX.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLEXPORT :: /PARAMS/, /DYNAUX/, /GRAPH/
#endif
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
      double precision NUH
      integer IEIGS
      double precision CONDC
      logical LSLACKS
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
      double precision pynrm, zpznrm, dnrm
      integer i
      character cslacks
      character*150 line(4)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( LSLACKS ) then
         cslacks = 's'
      else
         cslacks = ' '
      endif
      if( QFULL.eq.0 ) then
         pynrm = DNRM2(N, YPY, 1)
         zpznrm = DDOT(NIND, PZ, 1, PZ, 1)
C      if( QFULL.ne.1 ) then
C         zpznrm = zpznrm + DDOT(M, ZPZ, 1, ZPZ, 1)
C      endif
         if( zpznrm.gt.0.d0 ) zpznrm = DSQRT(zpznrm)
      else
         pynrm = 0.d0
         zpznrm = 0.d0
      endif
      dnrm = DNRM2(N, DX, 1)

      write(line,101)
      if( mod(ITER,10).eq.0 ) then
         call C_OUT(2,0,2,line)
      else
         call C_OUT(1,2,2,line)
      endif

      write(line,111) 
     1     ITER, ERR, C_ERR, MU, CNRM0, pynrm, cslacks, zpznrm,
     1     C_SKIP, dnrm, ALFAV, C_ALPHA, ALFAX, C_WATCH, NUH,
     1     LS_COUNT, F, IEIGS, CONDC
      call C_OUT(2,0,1,line)

C 101  format(/,'ITER     ERR        MU         C        ',
C     1           'YPY        ZPZ         D        ALFA1 ',
C     2           '     ALPHA  ',
C     2           '    NU     #LS       F       #eig  COND(C)')

 101  format(/,'ITER     ERR        MU         C        ',
     1           'YPY         PZ         D       ALFA(V)',
     2           '   ALFA(X)  ',
     2           '    NU     #LS       F       #cor  COND(C)')

 111  format(i5,d10.3,a1,d10.3,d10.3,d10.3,a1,d10.3,a1,d10.3,d10.3,a1,
     1       d10.3,a1,d10.3,i3,d16.8,i3,d11.3)
C 111  format(i5,G10.3,a1,G10.3,G10.3,G10.3,a1,G10.3,a1,G10.3,G10.3,a1,
C     1       G10.3,a1,G10.3,i3,G11.3,i3,G11.3)

C
C     If wished write even more info into a file
C
      if( QCNR.gt.0 .and. QPRINT.gt.1 ) then

         write(line,700) ITER
 700     format(/,80('-'),/,20(' '),
     1          'Detailed Information for Iteration ',
     1          i5,/,80('-'))
         call C_OUT(1,2,4,line)

         write(line,705) ITER,F/QFSCALE
 705     format(/,'  Objective function (Iter ',i5,'):',//,
     1        ' F = ',d25.17)
         call C_OUT(1,2,4,line)

         if( QPRINT.ge.9 ) then
            write(line,710) ITER
 710        format(/,'  Reduced bound multipliers (Iter ',i5,'):',/)
            call C_OUT(1,7,3,line)
            do i = 1, NIND
               write(line,711) i,RV(i)
 711           format(' RV(',i5,') = ',d25.17)
               call C_OUT(1,7,1,line)
            enddo
         endif

         if( QPRINT.ge.8 )then
            write(line,720) ITER
 720        format(/,'  Reduced gradient of objective function ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,6,3,line)
            do i = 1, NIND
               write(line,721) i,RG(i)
 721           format(' RG(',i5,') = ',d25.17)
               call C_OUT(1,6,1,line)
            enddo
         endif

         if( QPRINT.ge.8 )then
            write(line,730) ITER
 730        format(/,'  Reduced gradient of barrier function ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,6,3,line)
            do i = 1, NIND
               write(line,731) i,RGB(i)
 731           format(' RGB(',i5,') = ',d25.17)
               call C_OUT(1,6,1,line)
            enddo
         endif

         if( QPRINT.ge.7 .and. ITER.gt.0 )then
            write(line,740) ITER
 740        format(/,'  Quasi Newton Matrix ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,5,3,line)
            do i = 1, (NIND*(NIND+1))/2
               write(line,741) i,B(i)
 741           format(' B(',i5,') = ',d25.17)
               call C_OUT(1,5,1,line)
            enddo
         endif

         if( QPRINT.ge.8 )then
            write(line,750) ITER
 750        format(/,'  Sigmas for lower and upper bounds ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,6,3,line)
            do i = 1, NLB
               write(line,751) IVAR(ILB(i)),SIGMA_L(i)
 751           format(' SIGMA_L(',i5,') = ',d25.17)
               call C_OUT(1,6,1,line)
            enddo
            call C_OUT(1,6,0,line)
            do i = 1, NUB
               write(line,752) IVAR(IUB(i)),SIGMA_U(i)
 752           format(' SIGMA_U(',i5,') = ',d25.17)
               call C_OUT(1,6,1,line)
            enddo
         endif

         if( QPRINT.ge.7 .and. ITER.gt.0 )then
            write(line,760) ITER
 760        format(/,'  reduced Hessian ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,5,3,line)
            do i = 1, (NIND*(NIND+1))/2
               write(line,761) i,W(i)
 761           format(' W(',i5,') = ',d25.17)
               call C_OUT(1,5,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 ) then
            write(line,770) ITER
 770        format(/,'  range space step ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,4,3,line)
            do i = 1, N
               write(line,771) i,YPY(i)
 771           format(' YPY(',i5,') = ',d25.17)
               call C_OUT(1,4,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 )then
            write(line,780) ITER
 780        format(/,'  null space step (indepentend vars) ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,4,3,line)
            do i = 1, NIND
               write(line,781) i,PZ(i)
 781           format(' PZ(',i5,') = ',d25.17)
               call C_OUT(1,4,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 ) then
            write(line,790) ITER
 790        format(/,'  null space step (depentend vars) ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,4,3,line)
            do i = 1, M
               write(line,791) i,ZPZ(i)
 791           format(' ZPZ(',i5,') = ',d25.17)
               call C_OUT(1,4,1,line)
            enddo
         endif

         if( QPRINT.ge.6 .and. ITER.gt.0 )then
            write(line,795) ITER
 795        format(/,'  overall step ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,4,3,line)
            do i = 1, N
               write(line,796) IVAR(i),DX(i)
 796           format(' DX(',i5,') = ',d25.17)
               call C_OUT(1,4,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,800) ITER
 800        format(/,'  new primal variables ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,2,3,line)
            do i = 1, N
               write(line,801) IVAR(i),X(i)
 801           format(' X(',i5,') = ',d25.17)
               call C_OUT(1,2,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,810) ITER
 810        format(/,'  new slack variables ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,3,3,line)
            do i = 1, NLB
               write(line,811) IVAR(ILB(i)),S_L(i)
 811           format(' S_L(',i5,') = ',d25.17)
               call C_OUT(1,3,1,line)
            enddo
            call C_OUT(1,3,0,line)
            do i = 1, NUB
               write(line,812) IVAR(IUB(i)),S_U(i)
 812           format(' S_U(',i5,') = ',d25.17)
               call C_OUT(1,3,1,line)
            enddo
         endif

         if( QPRINT.ge.5 )then
            write(line,820) ITER
 820        format(/,'  new dual variables ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,3,3,line)
            do i = 1, NLB
               write(line,821) IVAR(ILB(i)),V_L(i)
 821           format(' V_L(',i5,') = ',d25.17)
               call C_OUT(1,3,1,line)
            enddo
            call C_OUT(1,3,0,line)
            do i = 1, NUB
               write(line,822) IVAR(IUB(i)),V_U(i)
 822           format(' V_U(',i5,') = ',d25.17)
               call C_OUT(1,3,1,line)
            enddo
         endif

         if( QPRINT.ge.5 .and. QLAMBDA.ne.0 )then
            write(line,825) ITER
 825        format(/,'  new multipliers for equalities ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,3,3,line)
            do i = 1, M
               write(line,826) i,LAM(i)
 826           format(' LAM(',i5,') = ',d25.17)
               call C_OUT(1,3,1,line)
            enddo
         endif

         if( QPRINT.ge.7 )then
            write(line,830) ITER
 830        format(/,'  complementarity: ',
     1             '(Iter ',i5,'):',/)
            call C_OUT(1,5,3,line)
            do i = 1, NLB
               write(line,831) IVAR(ILB(i)),IVAR(ILB(i)),V_L(i)*S_L(i)
 831           format(' V_L(',i5,')*S_L(',i5,') = ',d25.17)
               call C_OUT(1,5,1,line)
            enddo
            call C_OUT(1,5,0,line)
            do i = 1, NUB
               write(line,832) IVAR(IUB(i)),IVAR(IUB(i)),V_U(i)*S_U(i)
 832           format(' V_U(',i5,')*S_U(',i5,') = ',d25.17)
               call C_OUT(1,5,1,line)
            enddo
         endif

      endif
#ifdef WIN32
!--------------------------------------------------------------------
!      this section is for display iterative information graphically
!             YDL    03-06-2001

      if( OPT_FLAG ) then
         call my_draw_iter_info(%VAL(DRAWSUBP), iter,XORIG,ERR_RDGRD,
	1                 ERR_CNSTR,alfaV,alfaTau,alfaX,ERR,F, mu)
	endif
!--------------------------------------------------------------------
#endif
      return
      end

#ifdef WIN32
      subroutine my_draw_iter_info(DRAWSUB, iter,X,ERR_RDGRD,ERR_CNSTR,
	1                       alfaV,alfaTau,alfaX,ERR,F, mu)
C
C     Dirty trick to force draw_iter_info to be called from main executable
C
C     Andreas Waechter       10-04-01
C
      external DRAWSUB
	integer iter
	double precision x(*), err_rdgrd, err_cnstr, alfaV, alfaTau
      double precision alfaX, ERR, F, mu

	call DRAWSUB(iter,X,ERR_RDGRD,ERR_CNSTR,
	1             alfaV,alfaTau,alfaX,ERR,F, mu)
	return
	end
#endif
