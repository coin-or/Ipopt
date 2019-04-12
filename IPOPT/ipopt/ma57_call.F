C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA57_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS, FACT,
     1     LFACT, L_INC_LFACT, IFACT, LIFACT, L_INC_LIFACT, NEIGS,
     1     NZERO, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma57_call.F 615 2004-05-21 19:03:00Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA57
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
CP                          =1: perform symbolic manipulations with MA57AD
CP                          =2: factorize matrix from MA57BD
CP                          =3: use old factorization in order to solve for
CP                              rhs (MA57CD)
CP   NTOT      I    INT    total number of variables
CP   NZKKT     I    INT    number of nonzero elements in A
CP   A         I    DP     elements of matrix (nonly needed for TASK=2)
CP   IRN       I    INT    TASK=1: row indices of matrix
CP   JCN       I    INT    TASK=1: col indices of matrix
CP   RHS      I/O   DP     TASK=3: I: right hand sides
CP                                 O: solution
CP   FACT     I/O   DP     TASK=2: O:factors (need to be kept for TASK=3)
CP                         TASK=3: I:factors (from last TASK=2 call)
CP   LFACT    I/O   INT    TASK=1: O: Size of FACT need for later calls
CP                         TASK=2,3: I: Size of FACT provided
CP                         TASK=2: O: If too small (IERR=-1), new estimate
CP   L_INC_LFACT O  LOG    TASK=2: .true., if LFACT should be increased for
CP                                 next iteration
CP   IFACT    I/O   INT    TASK=2: O: IFACT for next call with TASK=3
CP                         TASK=3: I: IFACT from last call with TASK=2
CP   LIFACT   I/O   INT    TASK=1: O: Size of IFACT need for later calls
CP                         TASK=2,3: I: Size of IFACT provided
CP                         TASK=2: O: If too small (IERR=-1), new estimate
CP   L_INC_LIFACT O LOG    TASK=2: .true., if LIFACT should be increased for
CP                                 next iteration
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
CS    MA57ID
CS    MA57AD
CS    MA57BD
CS    MA57CD
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
      double precision A(NZKKT)
      integer IRN(NZKKT)
      integer JCN(NZKKT)
      double precision RHS(NTOT,*)
      integer LFACT
      double precision FACT(LFACT)
      logical L_INC_LFACT
      integer LIFACT
      integer IFACT(LIFACT)
      logical L_INC_LIFACT
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
      integer p_iwend, p_rwend, p_iwork
      integer i, iflag, job
      integer icntl(20), info(40)
      double precision cntl(5), rinfo(20)
      integer ierror, ncmpbr, ncmpbi

      integer P_KEEP, LKEEP
      save    P_KEEP, LKEEP

      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA57
      IERR = 97
      call C_OUT(2,0,1,
     1 'ma57_call:  Compiled without MA57 - cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         LKEEP  = 7*NTOT + NZKKT + MAX(NTOT,NZKKT) + 42
         P_KEEP = LIS
         LIS    = P_KEEP + LKEEP
         goto 9999
      endif
C
C     Normal execution
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Set parameters for MA57
C
      call MA57ID(cntl, icntl)
      icntl(1) = QCNR
      icntl(2) = QCNR
      icntl(3) = QCNR
      cntl(1)  = QPIVTOL
C
C     Jump to TASK
C
      goto (1000, 2000, 3000) TASK
      call C_OUT(2,0,1,'ma57_call: Wrong TASK!')
      stop
C
C     Get structure from MA57AD
C
 1000 continue
      p_iwork = p_iwend
      p_iwend = p_iwork + 5*NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      iflag = 0
      call MA57AD(NTOT, NZKKT, IRN, JCN, LKEEP, IS(P_KEEP+1),
     1     IW(p_iwork+1), icntl, info, rinfo)
      iflag  = info(1)
      ierror = info(2)
      if( iflag.ne.0 ) then
         write(line,*)
     1        'ma57_call: MA57AD returns info(1)=',info(1)
         call C_OUT(2,0,1,line)
         write(line,*)
     1        '           and info(2)=',info(2)
         call C_OUT(2,0,1,line)
         IERR = 508
         goto 9999
      endif
      p_iwend = p_iwork
C
      LFACT  = max(LFACT ,int(QFILLINFACT*dble(info(9))))
      LIFACT = max(LIFACT,int(QFILLINFACT*dble(info(10))))
C
      goto 9999
C
 2000 continue
C
C     Do the factorization
C
      p_iwork = p_iwend
      p_iwend = p_iwork + NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      call MA57BD(NTOT, NZKKT, A, FACT, LFACT, IFACT, LIFACT, LKEEP,
     1     IS(P_KEEP+1), IW(p_iwork+1), icntl, cntl, info, rinfo)
      iflag  = info(1)
      ierror = info(2)
      ncmpbr = info(28)
      ncmpbi = info(29)
      NEIGS  = info(24)
      if( iflag.eq.-4 ) then ! LIFACT too small
CTODO use fast update of MA57
         i = LIFACT
C         LIFACT = 10*info(18)
         LIFACT = int(QFILLINFACT*dble(info(18)))
         write(line,*) 'ma57_call:  LIFACT increased from',i,
     1        ' to',LIFACT
         call C_OUT(2,0,1,line)
         IERR = -1
         goto 9999
      elseif( iflag.eq.-3 ) then ! LFACT too small
CTODO use fast update of MA57
         i = LFACT
C         LFACT = 10*info(17)
         LFACT = int(QFILLINFACT*dble(info(17)))
         write(line,*) 'ma57_call:  LFACT increased from',i,
     1        ' to',LFACT
         call C_OUT(2,0,1,line)
         IERR = -1
         goto 9999
      elseif( iflag.eq.4 ) then   ! rank deficient matrix
         NZERO = 1
         NEIGS = 0
         goto 9999
      elseif( iflag.ne.0 ) then
         write(line,*)
     1        'ma57_call: MA57BD returns info(1) =',info(1)
         call C_OUT(2,0,1,line)
         write(line,*)
     1        '           and info(2)=',info(2)
         call C_OUT(2,0,1,line)
         IERR = 509
         goto 9999
      endif
C
C     Check if too many compressions were done; if so increase work space
C     for next iteration
C
      if( ncmpbr.ge.10 ) then
         write(line,*) 'ma57_call: NCMPBR=',ncmpbr,'; increase LFACT'
         call C_OUT(2,1,1,line)
         L_INC_LFACT = .true.
      else
         L_INC_LFACT = .false.
      endif
      if( ncmpbi.ge.10 ) then
         write(line,*) 'ma57_call: NCMPBI=',ncmpbi,
     1        '; increase LIFACT'
         call C_OUT(2,1,1,line)
         L_INC_LIFACT = .true.
      else
         L_INC_LIFACT = .false.
      endif
C
C     Determine inertia
C
      NZERO = 0
C
      goto 9999
C
C     Solve the system
C
 3000 continue
      p_iwork = p_iwend
      p_iwend = p_iwork + NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
      job = 1
CTODO: Check if multiple RHS Ok for MA57!
      call MA57CD(job, NTOT, FACT, LFACT, IFACT, LIFACT, NEIGS, RHS,
     1     NTOT, RW(p_rwend+1), LRW-p_rwend, IW(p_iwork+1), icntl, info)
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

      subroutine MA57_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
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
         LIW = 5*(N+M)
         LRW = 0
      elseif( TASK.eq.2 ) then
         LIW = (N+M)
         LRW = 0
      elseif( TASK.eq.3 ) then
         LIW = N+M
         LRW = 0
      endif


CAW   CHECK DOCUMENTATION ON LRW FOR MA57CD!
      LRW = N+M
      LIW = 5*(N+M)

      return
      end
