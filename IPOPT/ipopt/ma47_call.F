C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA47_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NEIGS, NZERO,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma47_call.F 632 2004-07-20 23:52:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA47
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
CP   TASK      I    INT    information about what to do:
CP                          =0: initialize, tell LRS_END, LIS_END
CP                          =1: Obtain structure from MA47
CP                          =2: factorize matrix by MA47BD
CP                          =3: use old factorization in order to solve for
CP                              rhs (MA47CD)
CP   NTOT      I    INT    total number of variables
CP   NZKKT     I    INT    number of nonzero elements in A
CP   A        I/O   DP     TASK=2: I:elements of matrix
CP                                 O:factors (need to be kept for TASK=3)
CP                         TASK=3: I:factors (from last TASK=2 call)
CP   IRN       I    INT    TASK=1: row indices of matrix
CP   JCN       I    INT    TASK=1: col indices of matrix
CP   RHS      I/O   DP     TASK=3: I: right hand sides
CP                                 O: solution
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
CS    MA47AD
CS    MA47BD
CS    MA47CD
CS    MA47ID
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
      integer p_iwend, p_rwend, p_w, p_iw1, p_irn
      integer icntl(7), info(24), i
      double precision cntl(2), rinfo(4)

      integer P_JCN, P_KEEP
      save    P_JCN, P_KEEP

      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA47
      IERR = 97
      call C_OUT(2,0,1,
     1 'ma47_call:  Compiled without MA47 - cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         P_JCN  = LIS
         P_KEEP = P_JCN + NZKKT
         LIS    = P_KEEP + NZKKT + 5*NTOT + 2
         goto 9999
      endif
C
C     Normal execution
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Get default values for MA47
C
      call MA47ID(cntl,icntl)
      icntl(1) = QCNR   ! redirect output of MA47 into file
      icntl(2) = QCNR   ! redirect output of MA47 into file
C      icntl(3) = 4      ! dump all output to file

      cntl(1) = QPIVTOL
C      cntl(2) = 1d-8
C
C     Jump to TASK
C
      goto (1000, 2000, 3000) TASK
      call C_OUT(2,0,1,'ma47_call: Wrong TASK!')
      stop
C
C     Get structure from MA47AD
C
 1000 continue
CTODO This could be different if overall KKT matrix does not need to be stored
C     for consistency check in get_step_full!
C     Note: Both IRN and JCN are changed by MA47AD
      p_irn   = p_iwend
      p_iwend = p_iwend + NZKKT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      do i = 1, NZKKT
         IW(p_irn+i) = IRN(i)
      enddo
      do i = 1, NZKKT
         IS(P_JCN+i) = JCN(i)
      enddo
      call MA47AD(NTOT, NZKKT, IW(p_irn+1), IS(P_JCN+1),
     1     IW(p_iwend+1), LIW-p_iwend, IS(P_KEEP+1), icntl,
     2     rinfo, info)
      if( info(1).lt.0 .or. info(1).ge.4 ) then
         write(line,*)
     1        'ma47_call: MA47AD returns INFO(1)=',INFO(1)
         call C_OUT(2,0,1,line)
         IERR = 510
         goto 9999
      endif
      p_iwend = p_irn
C
      LIWMA = max(LIWMA,int(QFILLINFACT*dble(info(7))))
      LA    = max(LA,int(QFILLINFACT*dble(info(6))))
C
      goto 9999
C
 2000 continue
C
C     Do the factorization
C
      p_iw1   = p_iwend
      p_iwend = p_iw1 + 2*NTOT + 2
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      call MA47BD(NTOT, NZKKT, IS(P_JCN+1), A, LA, IWMA,
     1     LIWMA, IS(P_KEEP+1), cntl, icntl, IW(p_iw1+1), rinfo, info)
      if( info(1).lt.0 .or. info(1).ge.4 ) then
         if( info(1).eq.-3 ) then   ! LIWMA too small
            i = LIWMA
C            LIWMA = 10*info(2)
            LIWMA = int(QFILLINFACT*dble(info(2)))
            write(line,*) 'ma47_call:  LIWMA increased from',i,
     1           ' to',LIWMA
            call C_OUT(2,1,1,line)
            IERR = -1
            goto 9999
         endif
         if( info(1).eq.-4 ) then   ! LA too small
            i = LA
C            LA = 10*info(2)
            LA = int(QFILLINFACT*dble(info(2)))
            write(line,*) 'ma47_call:  LA increased from',i,
     1           ' to',LA
            call C_OUT(2,1,1,line)
            IERR = -1
            goto 9999
         endif
         write(line,*)
     1        'ma47_call: MA47BD returns INFO(1)=',INFO(1)
         call C_OUT(2,0,1,line)
         IERR = 511
         goto 9999
      endif
C
C     Check if too many compressions were done; if so increase work space
C     for next iteration
C
      if( INFO(18).ge.10 ) then
         write(line,*) 'ma47_call: INFO(18)=',INFO(18),'; increase LA'
         call C_OUT(2,1,1,line)
         L_INC_LA = .true.
      else
         L_INC_LA = .false.
      endif
      if( INFO(19).ge.10 ) then
         write(line,*) 'ma47_call: INFO(19)=',INFO(19),
     1        '; increase LIWMA'
         call C_OUT(2,1,1,line)
         L_INC_LIWMA = .true.
      else
         L_INC_LIWMA = .false.
      endif
C
C     Determine inertia
C
      NEIGS = INFO(23)
      NZERO = INFO(24)
C
      goto 9999
C
C     Solve the system
C
 3000 continue
      p_w     = p_rwend
      p_rwend = p_w + NTOT
      p_iw1   = p_iwend
      p_iwend = p_iw1 + NTOT
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
      do i = 1, NEIGS
         call MA47CD(NTOT, A, LA, IWMA, LIWMA, RW(p_w+1),
     1        RHS(1,i), IW(p_iw1+1), icntl)
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

      subroutine MA47_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
     1     LRW, LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer TASK, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT

      if( TASK.eq.0 ) then
         ! Estimate LA and LIWMA
         LIW = int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
         LRW = int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
      elseif( TASK.eq.1 ) then
         LIW = NZKKT + ((4+2*NZKKT+5*(N+M))*3)/2
         LRW = 0
      elseif( TASK.eq.2 ) then
         LIW = 2*(N+M) + 2
         LRW = 0
      elseif( TASK.eq.3 ) then
         LIW = N+M
         LRW = N+M
      endif

      return
      end
