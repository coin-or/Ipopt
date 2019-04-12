C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine WGSMP_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NEIGS, NZERO,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: wgsmp_call.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to WGSMP
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
CA    Andreas Waechter      02/16/03  Initial version
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
CS    HSL2GCSC
CS    WGSMP
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
      integer p_iwend, p_rwend, p_avals, i, j
      double precision dummy

      integer IPARM(64)
      double precision DPARM(64)
      save    IPARM,   DPARM

      integer P_IA, P_JA, P_PHSL2GCSC, P_PHSL2GCSC2, NZGCSC
      save    P_IA, P_JA, P_PHSL2GCSC, P_PHSL2GCSC2, NZGCSC

      logical HAVE_HSL2GCSC, HAVE_ANALYSE
      save    HAVE_HSL2GCSC, HAVE_ANALYSE

      character*80 line(2)

      integer dump, wantcount
      parameter( dump = 0, wantcount = 1 )
      integer count
      save count
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
CTODO CHECK ALL!
 9901 format(i16)
 9910 format(d23.15)
#ifndef HAVE_WSMP
      IERR = 97
      call C_OUT(2,0,1,
     1'wgsmp_call:  Compiled without WSMP - cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         P_IA         = LIS
         P_JA         = P_IA         + NTOT + 1
         P_PHSL2GCSC  = P_JA         + 2*NZKKT
         P_PHSL2GCSC2 = P_PHSL2GCSC  + NZKKT
         LIS          = P_PHSL2GCSC2 + NZKKT

         if( dump.ne.0 ) then
            write(45,9901) NZKKT
            count = 0
         endif
         HAVE_HSL2GCSC = .false.

         goto 9999
      endif
C
C     Normal execution
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0
C
C     Jump to TASK
C
      goto (1000, 2000, 3000) TASK
      call C_OUT(2,0,1,'wgsmp_call: Wrong TASK!')
      stop
C
C     Get structure from MA27AD
C
 1000 continue
C
C     First check if we have to get the HSL format to Compressed Sparse
C     Column (CSC) format
C
      if( .not.HAVE_HSL2GCSC ) then

         HAVE_HSL2GCSC = .true.
         NZGCSC = 2*NZKKT
         call HSL2GCSC(NTOT, NZKKT, IRN, JCN, NZGCSC, IS(P_IA+1),
     1        IS(P_JA+1), IS(P_PHSL2GCSC+1), IS(P_PHSL2GCSC2+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*) 'wsmp_call: HSL2GCSC returned IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 133
            goto 9999
         endif
C
C     Get default values for WSMP
C
         IPARM(1) = 0
         IPARM(2) = 0
         IPARM(3) = 0
         call WGSMP( NTOT, IS(P_IA+1), IS(P_JA+1), dummy, dummy,
     1        0, 0, dummy, IPARM, DPARM)
C
C     Set them to what we need lateron
C
c          IPARM(8) = 1
c          DPARM(11) = 1.d-8  ! Pivot tolernace
          IPARM(16) = 1 ! Take long time for ordering phase
c          IPARM(17) = NTOT+1 ! pure minimum local fill-in
c          IPARM(18) = 1      ! even less loval fill in


C         IPARM(15) = 0 ! no restrictions on pivoting
C         IPARM(31) = 1 ! want L D L^T decomposition
C         IPARM(32) = 1 ! need D to see where first negative eigenvalue occurs
C
C         IPARM(20) = 1 ! for ordering in IP methods?
C
C         IPARM(11) = 1
C         DPARM(11) = 4.d-16
C         DPARM(21) = 1.d+200

       endif

      HAVE_ANALYSE = .false.

      LA = NZKKT
      LIWMA = 0

      goto 9999

 2000 continue
C
C     Do the factorization
C
      count = count + 1
      if( HAVE_ANALYSE ) then
         IPARM(2) = 2
      else
         IPARM(2) = 1
         HAVE_ANALYSE = .true.
      endif
      IPARM(3) = 2

      p_avals = p_rwend
      p_rwend = p_avals + NZGCSC
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Copy elements into CSC format
C
      call DCOPY(NZGCSC, 0.D0, 0, RW(p_avals+1), 1)
      do i = 1, NZKKT
         j = IS(P_PHSL2GCSC+i)
         RW(p_avals+j) = RW(p_avals+j) + A(i)
      enddo
      do i = 1, NZKKT
         j = IS(P_PHSL2GCSC2+i)
         if( j.ne.-1 ) RW(p_avals+j) = RW(p_avals+j) + A(i)
      enddo
      if( dump.eq.2 .and. count.eq.wantcount ) then
         do i = 1, NTOT
            do j = IS(P_IA+i), IS(P_IA+i+1)-1
               write(44,*) 'A(',i,',',IS(P_JA+j),')=',RW(p_avals+j)
            enddo
         enddo
      endif
      if( dump.eq.2 ) then
         do i = 1, 64
            write(45,9901) IPARM(i)
         enddo
         do i = 1, 64
            write(45,9910) DPARM(i)
         enddo
      endif
      if( dump.eq.2 .and. count.eq.wantcount ) then
         write(45,*) 'AVALS'
         do i = 1, IS(P_IA+NTOT+1)-1
            write(45,9910) RW(p_avals+i)
         enddo
         write(45,*) 'AVALS End'
      endif
      if( dump.eq.2 .and. count.eq.wantcount ) then
         write(*,*) 'dump = ',dump
         stop
      endif
      NZERO = 0
      call WGSMP(NTOT, IS(P_IA+1), IS(P_JA+1), RW(p_avals+1),
     1     dummy, 0, 0, dummy, IPARM, DPARM)
      if( IPARM(64).gt.0 ) then
         NZERO = 1
      elseif( IPARM(64).ne.0 ) then
         write(line,*)
     1        'wgsmp_call:  Error accurred in WGSMP-2 IPARM(64)=',
     2        IPARM(64)
         call C_OUT(2,0,1,line)
         if( IPARM(64).eq.-102 ) then
            IERR = 96
         else
            IERR = 592
         endif
         goto 9900
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'wgsmp_call:  Size of work space (analysis) = ',
     1        IPARM(24)
         call C_OUT(2,0,1,line)
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'wgsmp_call:  Size of work space (factor) = ',
     1        IPARM(23)
         call C_OUT(2,0,1,line)
      endif

      NEIGS = 0
      L_INC_LA    = .false.
      L_INC_LIWMA = .false.

      goto 9999
C
C     Solve the system
C
 3000 continue
      IPARM(2) = 3
      IPARM(3) = 3              ! No iterative refinement
      IPARM(3) = 4              ! With iterative refinement
      if( IPARM(3).eq.4 ) then
         p_avals = p_rwend
         p_rwend = p_avals + NZGCSC
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     Copy elements into CSC format
C
         call DCOPY(NZGCSC, 0.D0, 0, RW(p_avals+1), 1)
         do i = 1, NZKKT
            j = IS(P_PHSL2GCSC+i)
            RW(p_avals+j) = RW(p_avals+j) + A(i)
         enddo
         do i = 1, NZKKT
            j = IS(P_PHSL2GCSC2+i)
            if( j.ne.-1 ) RW(p_avals+j) = RW(p_avals+j) + A(i)
         enddo
      else
         p_avals = 0
      endif
      call WGSMP(NTOT, IS(P_IA+1), IS(P_JA+1), RW(p_avals+1), RHS, NTOT,
     1     NEIGS, dummy, IPARM, DPARM)
      if( IPARM(64).ne.0 ) then
         write(line,*)
     1        'wgsmp_call:  Error accurred in WSMP-3-4 IPARM(64)=',
     2        IPARM(64)
         call C_OUT(2,0,1,line)
         if( IPARM(64).eq.-102 ) then
            IERR = 96
         else
            IERR = 592
         endif
         goto 9900
      endif

      goto 9999
C
 9900 continue
      if( QPRINT.ge.3 ) then
         do i = 1, 64
            write(line,*) 'IPARM(',i,') = ',IPARM(i)
            call C_OUT(1,1,1,line)
         enddo
         do i = 1, 64
            write(line,*) 'DPARM(',i,') = ',DPARM(i)
            call C_OUT(1,1,1,line)
         enddo
      endif
#endif
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine WGSMP_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
     1     LRW, LIW, FILLINFACT)

      implicit none
      include 'IPOPT.INC'
      integer TASK, N, M, NLB, NUB, NZA, NZKKT, LRW, LIW
      double precision FILLINFACT

      if( TASK.eq.0 ) then
         ! Estimate LA and LIWMA
         LIW = 0
         LRW = 0
      elseif( TASK.eq.1 ) then
         LIW = 0
         LRW = 2*NZKKT
      elseif( TASK.eq.2 ) then
         LIW = 0
         LRW = 2*NZKKT
      elseif( TASK.eq.3 ) then
         LIW = 0
         LRW = 2*NZKKT
      endif

      return
      end
