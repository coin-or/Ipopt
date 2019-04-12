C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine MA27_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NEIGS, NZERO,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: ma27_call.F 675 2005-07-26 18:47:13Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to MA27
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
CP                          =1: Obtain structure from MA27
CP                          =2: factorize matrix by MA27BD
CP                          =3: use old factorization in order to solve for
CP                              rhs (MA27CD)
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
CS    MA27AD
CS    MA27BD
CS    MA27CD
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
C*******************************************************************************
C
C                              COMMON blocks
C
C*******************************************************************************
C
#ifdef HAVE_MA27
#ifdef OLD_MA27
      double precision U
      integer LP, MP, LDIAG
      common /MA27DD/ U, LP, MP, LDIAG

      double precision OPS
      integer IERROR, NRLTOT, NIRTOT, NRLNEC, NIRNEC, NRLADU, NIRADU
      integer NRLBDU, NIRBDU, NCMPA, NCMPBR, NCMPBI, NTWO, NEIG
      common /MA27ED/ OPS, IERROR, NRLTOT, NIRTOT, NRLNEC, NIRNEC,
     1     NRLADU, NIRADU, NRLBDU, NIRBDU, NCMPA, NCMPBR, NCMPBI,
     2     NTWO, NEIG
#endif
#endif
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_iwend, p_rwend, p_w, p_iw1
      integer i, iflag
#ifndef OLD_MA27
      integer icntl(30), info(20)
      double precision cntl(5), ops
      integer ierror, nirnec, nrlnec, ncmpbr, ncmpbi, neig
#endif

      integer P_IKEEP, NSTEPS, MAXFRT
      save    P_IKEEP, NSTEPS, MAXFRT

      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
#ifndef HAVE_MA27
      IERR = 97
      call C_OUT(2,0,1,
     1 'ma27_call:  Compiled without MA27 - cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         P_IKEEP = LIS
         LIS     = P_IKEEP + 3*NTOT
         goto 9999
      endif
C
C     Normal execution
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Set default values for MA27
C
#ifdef OLD_MA27
      LP = QCNR
      MP = QCNR
      U = QPIVTOL
#else
      call MA27ID(icntl, cntl)
      icntl(1) = QCNR
      icntl(2) = QCNR
      cntl(1)  = QPIVTOL
#endif
C
C     Jump to TASK
C
      goto (1000, 2000, 3000) TASK
      call C_OUT(2,0,1,'ma27_call: Wrong TASK!')
      stop
C
C     Get structure from MA27AD
C
 1000 continue
CTODO This could be different if overall KKT matrix does not need to be stored
C     for consistency check in get_step_full!
C     Note: Both IRN and JCN are changed by MA27AD
      p_iw1   = p_iwend
      p_iwend = p_iw1 + 2*NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      iflag = 0
#ifdef OLD_MA27
      call MA27AD(NTOT, NZKKT, IRN, JCN, IW(p_iwend+1), LIW-p_iwend,
     1     IS(P_IKEEP+1), IW(p_iw1+1), NSTEPS, iflag)
#else
      call MA27AD(NTOT, NZKKT, IRN, JCN, IW(p_iwend+1), LIW-p_iwend,
     1     IS(P_IKEEP+1), IW(p_iw1+1), NSTEPS, iflag, icntl, cntl,
     2        info, ops)
      iflag  = info(1)
      nrlnec = info(5)
      nirnec = info(6)
#endif
      if( iflag.ne.0 ) then
         write(line,*)
     1        'ma27_call: MA27AD returns iflag=',iflag
         call C_OUT(2,0,1,line)
         IERR = 508
         if( iflag.eq.-3 ) IERR = 99
         goto 9999
      endif
      p_iwend = p_iw1
C
      LIWMA = max(LIWMA,int(QFILLINFACT*dble(NIRNEC)))
      LA    = max(LA   ,int(QFILLINFACT*dble(NRLNEC)))
C
      goto 9999
C
 2000 continue
C
C     Do the factorization
C
      p_iw1   = p_iwend
      p_iwend = p_iw1 + NTOT
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
c$$$CDELETEME
c$$$      write(QCNR,*) 'MA27 matrix'
c$$$      do i = 1, NZKKT
c$$$         write(QCNR,10024) i, IRN(i), JCN(i), A(i)
c$$$10024    format(i6,i6,i6,d25.16)
c$$$      enddo
#ifdef OLD_MA27
      call MA27BD(NTOT, NZKKT, IRN, JCN, A, LA, IWMA, LIWMA,
     1     IS(P_IKEEP+1), NSTEPS, MAXFRT, IW(p_iw1+1), iflag)
#else
      call MA27BD(NTOT, NZKKT, IRN, JCN, A, LA, IWMA, LIWMA,
     1     IS(P_IKEEP+1), NSTEPS, MAXFRT, IW(p_iw1+1), icntl,
     2     cntl, info)
      iflag  = info(1)
      ierror = info(2)
      ncmpbr = info(12)
      ncmpbi = info(13)
      neig   = info(15)
#endif
      if( iflag.eq.-3 .or. iflag.eq.-4) then ! LIWMA or LA too small
         i = LIWMA
C         LIWMA = 10*IERROR
         LIWMA = int(QFILLINFACT*dble(IERROR))
         write(line,*) 'ma27_call:  LIWMA increased from',i,
     1        ' to',LIWMA
         call C_OUT(2,1,1,line)
         i = LA
C         LA = 10*IERROR
         LA = int(QFILLINFACT*dble(IERROR))
         write(line,*) 'ma27_call:  LA increased from',i,
     1        ' to',LA
         call C_OUT(2,1,1,line)
         IERR = -1
         goto 9999
      elseif( iflag.eq.-5 .or. iflag.eq.3 ) then   ! singular matrix
         NZERO = 1
         NEIGS = 0
         goto 9999
      elseif( iflag.ne.0 ) then
         write(line,*)
     1        'ma27_call: MA27BD returns iflag =',iflag
         call C_OUT(2,0,1,line)
         IERR = 509
         goto 9999
      endif
C
C     Check if too many compressions were done; if so increase work space
C     for next iteration
C
      if( NCMPBR.ge.10 ) then
         write(line,*) 'ma27_call: NCMPBR=',NCMPBR,'; increase LA'
         call C_OUT(2,1,1,line)
         L_INC_LA = .true.
      else
         L_INC_LA = .false.
      endif
      if( NCMPBI.ge.10 ) then
         write(line,*) 'ma27_call: NCMPBI=',NCMPBI,
     1        '; increase LIWMA'
         call C_OUT(2,1,1,line)
         L_INC_LIWMA = .true.
      else
         L_INC_LIWMA = .false.
      endif
C
C     Determine inertia
C
      NEIGS = NEIG
      NZERO = 0
C
      goto 9999
C
C     Solve the system
C
 3000 continue
      p_w     = p_rwend
      p_rwend = p_w + MAXFRT
      p_iw1   = p_iwend
      p_iwend = p_iw1 + NSTEPS
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
      do i = 1, NEIGS
#ifdef OLD_MA27
         call MA27CD(NTOT, A, LA, IWMA, LIWMA, RW(p_w+1), MAXFRT,
     1        RHS(1,i), IW(p_iw1+1), NSTEPS)
#else
         call MA27CD(NTOT, A, LA, IWMA, LIWMA, RW(p_w+1), MAXFRT,
     1        RHS(1,i), IW(p_iw1+1), NSTEPS, icntl, cntl)
#endif
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

      subroutine MA27_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
     1     LRW, LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer TASK, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT
      character*80 line

      if( TASK.eq.0 ) then
         ! Estimate LA and LIWMA
         LIW = int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
         LRW = int(QFILLINFACT*FILLINFACT*(dble(NZKKT)))
      elseif( TASK.eq.1 ) then
         LIW = 2*(N+M) + ((1+2*NZKKT+3*(N+M))*3)/2
         LRW = 0
      elseif( TASK.eq.2 ) then
         LIW = N+M
         LRW = 0
      elseif( TASK.eq.3 ) then
         LIW = N+M
         LRW = N+M
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'ma27_call_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
