C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA48SYM_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS,
     1     LIWMA, L_INC_LIWMA, IWMA, LA, L_INC_LA, NEIGS, NZERO,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma48sym_call.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA48 for solving SYMMETRIC systems
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
CA    Andreas Waechter      11/16/02   recovered from earlier archive
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
CP   TASK      I    INT    information about what to do:
CP                          =0: initialize, tell LRS_END, LIS_END
CP                          =1: Obtain structure from MA48
CP                          =2: factorize matrix by MA48BD
CP                          =3: use old factorization in order to solve for
CP                              rhs (MA48CD)
CP   NTOT      I    INT    total number of variables
CP   NZKKT     I    INT    number of nonzero elements in A
CP   A        I/O   DP     TASK=2: I:elements of matrix
CP                                 O:factors (need to be kept for TASK=3)
CP                         TASK=3: I:factors (from last TASK=2 call)
CP   IRN       I    INT    TASK=1: row indices of matrix
CP   JCN       I    INT    TASK=1: col indices of matrix
CP   RHS      I/O   DP     TASK=3: I: right hand sides
CP                                 O: solutions
CP   LIWMA    I/O   LOG    TASK=1: O: Size of IWMA need for later calls
CP                         TASK=2,3: I: Size of IWMA provided
CP                         TASK=2: O: If too small (IERR=-1), new estimate
CP   L_INC_LIWMA O  INT    TASK=2: .true., if LIWMA should be increased for next
CP                                 iteration
CP   IWMA     I/O   INT    TASK=2: O: IW for next call with TASK=3
CP                         TASK=3: I: IW from last call with TASK=2
CP   LA       I/O   INT    TASK=1: O: Size LA of A need for later calls
CP                         TASK=2,3: I: Size LA of A provided
CP                         TASK=2: O: If too small (IERR=-1), new estimate
CP   L_INC_LA  O    LOG    TASK=2: .true., if LA should be increased for next
CP                                 iteration
CP   NEIGS     O    INT    TASK=2: number of negative eigenvalues
CP             I    INT    TASK=3: number of right hand sides
CP   NZERO     O    INT    TASK=2: number of zero eigenvalues
CP   LRS      I/O   INT    length of RS (Output for TASK = 0)
CP   RS       I/O   DP     can be used to store DP variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LIS      I/O   INT    length of IS (Output for TASK = 0)
CP   IS       I/O   INT    can be used to store INT variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LIW       I    INT    length of IW (Output for TASK = 0)
CP   IW        W    INT    integer work space
CP   LRW       I    INT    length of RW (Output for TASK = 0)
CP   RW        W    DP     double precision work space
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         =-1: LIWMA or LA are too small, start over with
CP                              restored matrix
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
CS    MA48AD
CS    MA48BD
CS    MA48CD
CS    MA48ID
CS    C_OUT
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
      integer TASK
      integer NTOT
      integer NZKKT
      double precision A(*)
      integer IRN(*)
      integer JCN(*)
      double precision RHS(NTOT,*)
      integer LIWMA
      logical L_INC_LIWMA
      integer IWMA(LIWMA)
      integer LA
      logical L_INC_LA
      integer NEIGS
      integer NZERO
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_iwend, p_rwend, p_w, p_iw, p_x
      integer icntl(9), info(12), i, job
      double precision cntl(5), rinfo, error(3)
      logical transp_l

      integer P_IRN, P_JCN, P_KEEP, P_A, LA_INT, NZKKT_INT, JOB_B
      save    P_IRN, P_JCN, P_KEEP, P_A, LA_INT, NZKKT_INT, JOB_B

      logical COPIED_A
      save    COPIED_A

      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA48
      IERR = 97
      call C_OUT(2,0,1, 'ma48sym_call:  Compiled without MA48 - '
     1     //'cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         LA_INT = max(2*NZKKT,int(QFILLINFACT*NZKKT)) !This is to store factor INSIDE this routine
         P_IRN  = LIS
         P_JCN  = P_IRN + LA_INT
         P_KEEP = P_JCN + LA_INT
         LIS    = P_KEEP + 10*NTOT + 7
         P_A    = LRS
         LRS    = P_A + LA_INT
         COPIED_A = .false.
         JOB_B = 1
         goto 9999
      endif
C
C     Normal execution
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Get default values for MA48
C
      call MA48ID(cntl,icntl)
      icntl(1) = QCNR   ! redirect output of MA48 into file
      icntl(2) = QCNR   ! redirect output of MA48 into file
C      icntl(3) = 4      ! dump all output to file
      cntl(2) = QPIVTOL
C
C     Jump to TASK
C
      goto (1000, 2000, 3000) TASK
      call C_OUT(2,0,1,'ma48sym_call: Wrong TASK!')
      stop
C
C     Get structure from MA48AD
C
 1000 continue
CTODO
      call C_OUT(2,0,1,'WARNING: MA48 option not well tested!')

CTODO This could be different if overall KKT matrix does not need to be stored
C     for consistency check in get_step_full!
      NZKKT_INT = 0
      do i = 1, NZKKT
         if( IRN(i).ne.JCN(i) ) then
            IS(P_IRN+NZKKT_INT+1) = IRN(i)
            IS(P_JCN+NZKKT_INT+1) = JCN(i)
            RS(P_A  +NZKKT_INT+1) = A(i)
            IS(P_IRN+NZKKT_INT+2) = JCN(i)
            IS(P_JCN+NZKKT_INT+2) = IRN(i)
            RS(P_A  +NZKKT_INT+2) = A(i)
            NZKKT_INT = NZKKT_INT + 2
         else
            NZKKT_INT = NZKKT_INT + 1
            IS(P_IRN+NZKKT_INT) = IRN(i)
            IS(P_JCN+NZKKT_INT) = JCN(i)
            RS(P_A  +NZKKT_INT) = A(i)
         endif
      enddo
      p_iw    = p_iwend
      p_iwend = p_iw + 9*NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      if( QSYMSOLV.eq.5 ) then
         job = 1
      else
         job = 3
      endif
      call MA48AD(NTOT, NTOT, NZKKT_INT, job, LA_INT, RS(P_A+1),
     1     IS(P_IRN+1), IS(P_JCN+1), IS(P_KEEP+1), cntl, icntl,
     1     IW(p_iw+1), info, rinfo )
      if( info(1).lt.0 ) then
         write(line,213) info(1)
 213     format ('ma48sym_call: MA48AD returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-3 ) then
            call C_OUT(2,0,1,
     1           'ma48sym_call: Value of DFILLINFACT too small.')
         elseif( info(1).eq.-4 ) then
            write(line,214) info(10)
 214        format('   Matrix is structurally rank deficient with rank',
     1           i4)
            call C_OUT(2,0,1,line)
         endif
         IERR = 528
         goto 9999
      elseif( info(1).eq.2 ) then
         NZERO = NTOT - info(5)
         goto 9999
      elseif( info(1).gt.0 ) then
         write(line,211) info(1)
 211     format('ma48sym_call: MA48AD returns info(1)=',i2)
         call C_OUT(2,1,1,line)
      endif
      if( info(2).gt.10 ) then
         write(line,212) info(2)
 212     format('ma48sym_call: Number of garbage collections=',i2,/,
     1          '           LA should be increased!')
         call C_OUT(2,0,2,line)
      endif
C
C     This is a flag to indicate that for the next call of MA48BD we don't need
C     to copy A again...
C
      COPIED_A = .false.
C
C     Don't need extra large arrays here, since everything is stored internally
C     for now
C
      LIWMA = 0
      LA    = NZKKT
C
      goto 9999
C
 2000 continue
C
C     Get work space pointers
C
      p_w     = p_rwend
      p_rwend = p_w     + NTOT
      p_iw    = p_iwend
      p_iwend = p_iw    + 4*NTOT
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Don't need this when MA48AD has been called before...
C
      if( .not.COPIED_A ) then
         NZKKT_INT = 0
         do i = 1, NZKKT
            if( IRN(i).ne.JCN(i) ) then
               RS(P_A  +NZKKT_INT+1) = A(i)
               RS(P_A  +NZKKT_INT+2) = A(i)
               NZKKT_INT = NZKKT_INT + 2
            else
               NZKKT_INT = NZKKT_INT + 1
               RS(P_A  +NZKKT_INT) = A(i)
            endif
         enddo
      endif
      COPIED_A = .false.
C
C     Do the factorization
C
      call MA48BD(NTOT, NTOT, NZKKT_INT, JOB_B, LA_INT, RS(P_A+1),
     1     IS(P_IRN+1), IS(P_JCN+1), IS(P_KEEP+1), cntl, icntl,
     2     RW(p_w+1), IW(p_iw+1), info, rinfo)
      if( info(1).eq.2 ) then
         NZERO = 1
         goto 9999
      elseif( info(1).ne.0 ) then
C      if( info(1).lt.0 ) then
         write(line,251) info(1)
 251     format ('ma48sym_call: MA48BD returns info(1)=',i3)
         call C_OUT(2,0,1,line)
         if( info(1).eq.-3 ) then
            call C_OUT(2,0,1,
     1           'ma48sym_call: Value of DFILLINFACT too small.')
         endif
         IERR = 524
         goto 9999
      endif
CMOD: If you uncomment the next line, MA48BC will use the fast factorization
C     next time, but this might be less accurate and lead to more iterations
c      JOB_B = 2

C
C     Check if too many compressions were done; if so increase work space
C     for next iteration
C
      L_INC_LA = .false.
      L_INC_LIWMA = .false.
C
C     Determine inertia
C
      NEIGS = 0                 ! Change this to M in get_step_full
      NZERO = 0
C
      goto 9999
C
C     Solve the system
C
 3000 continue
C
C     Get work space pointers
C
      p_x     = p_rwend
      p_w     = p_x     + NTOT
      p_rwend = p_w     + 4*NTOT
      p_iw    = p_iwend
      p_iwend = p_iw    + NTOT
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      transp_l = .false.
CTODO The following choice determines whether MA48CD should use iterative
C     retinement.  If job = 1, then no iterative refinement is performed,
C     and the calling routine might perform iterative refinement directly.
C     If job = 3, then MA48CD performs iterative refinement, which seems to
C     perform quite a bit slower, and it also seemed to me that it sometime
C     gave bad solutions.
c      job = 3                   ! for now, with iterative refinement
      job = 1                   ! no iterative refinement
      do i = 1, NEIGS
         call MA48CD(NTOT, NTOT, transp_l, job, LA_INT, RS(P_A+1),
     1        IS(P_IRN+1), IS(P_KEEP+1), cntl, icntl, RHS(1,i),
     2        RW(p_x+1), error, RW(p_w+1), IW(p_iw+1), info)
         if( info(1).ne.0 ) then
            write(line,301) info(1)
 301        format ('ma48sym_call: MA48CD returns info(1)=',i3)
            call C_OUT(2,0,1,line)
            IERR = 526
            goto 9999
         endif
         call DCOPY(NTOT, RW(p_x+1), 1, RHS, 1)
      enddo
      goto 9999
C
#endif
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine MA48SYM_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
     1     LRW, LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer TASK, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT

      if( TASK.eq.0 ) then
         ! Estimate LA and LIWMA
         LIW = 2*int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
         LRW = 2*int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
      elseif( TASK.eq.1 ) then
         LIW = 9*(N+M)
         LRW = 0
      elseif( TASK.eq.2 ) then
         LIW = 4*(N+M)
         LRW = N+M
      elseif( TASK.eq.3 ) then
         LIW = N+M
         LRW = 5*(N+M)
      endif

      return
      end
