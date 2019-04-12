C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine WSSMP_CALL(TASK, NTOT, NZKKT, A, IRN, JCN, RHS, LIWMA,
     1     L_INC_LIWMA, IWMA, LA, L_INC_LA, NEIGS, NZERO,
     1     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: wssmp_call.F 615 2004-05-21 19:03:00Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface to WSSMP
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
CA    Andreas Waechter      06/11/02  Initial version
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
CP             I           TASK=1: number of variables (fix for WSSMP)
CP                                 (only needed if pivoting is restricted!)
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
CP                           On first call with TASK=1 and QSYMSOLV=4, the
CP                           first N (=neigs) entries are 1 if corresponding
CP                           diagonal element in upper right block of iteration
CP                           matrix not structually zero.
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
CS    HSL2CSC
CS    WSETMAXTHRDS
CS    WSSMP
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
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer p_iwend, p_rwend, p_diag, p_avals
      integer i, idummy, j
      double precision dummy

      integer IPARM(64)
      double precision DPARM(64)
      save    IPARM,   DPARM

      integer P_IA, P_JA, P_PHSL2CSC, P_PERM, P_INVP, NZCSC
      save    P_IA, P_JA, P_PHSL2CSC, P_PERM, P_INVP, NZCSC

      logical HAVE_HSL2CSC
      save    HAVE_HSL2CSC

      character*80 line(2)

C     the following can be used to dump the matrix to a file
      integer dump, wantcount
      parameter( dump = 0, wantcount = 10 )
      integer count, nthreads, nthreadsmax, nvarperthread
      parameter( nthreadsmax = 32, nvarperthread = 1000 )
      save count
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
 9901 format(2i8)
 9910 format(d23.16,i16)
#ifndef HAVE_WSMP
      IERR = 97
      call C_OUT(2,0,1,
     1'wssmp_call:  Compiled without WSMP - cannot execute this option')
      goto 9999
#else
C
C     Initialize storage space
C
      if( TASK.eq.0) then
         P_IA       = LIS
         P_JA       = P_IA       + NTOT + 1
         P_PHSL2CSC = P_JA       + NZKKT ! here we assume that
                                         ! all diagonal are included!
         P_PERM     = P_PHSL2CSC + NZKKT
         P_INVP     = P_PERM     + NTOT
         LIS        = P_INVP     + NTOT

         if( dump.ne.0 ) then
            write(45,9901) NZKKT
            count = 0
         endif
         HAVE_HSL2CSC = .false.

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
      call C_OUT(2,0,1,'wssmp_call: Wrong TASK!')
      stop
C
C     Get structure from MA27AD
C
 1000 continue
C
C     First check if we have to get the HSL format to Compressed Sparse
C     Column (CSC) format
C
      if( .not.HAVE_HSL2CSC ) then

         HAVE_HSL2CSC = .true.
         NZCSC = NZKKT
         call HSL2CSC(NTOT, NZKKT, IRN, JCN, NZCSC, IS(P_IA+1),
     1        IS(P_JA+1), IS(P_PHSL2CSC+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*) 'wssmp_call: HSL2CSC returned IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 133
            goto 9999
         endif
C
C     The followoing is a trick to speed up factorization is we have no
C     supernodes (more threads than processors...)
C
         nthreads = min( nthreadsmax, NTOT/nvarperthread + 1)
         if( QPRINT.ge.2 ) then
            write(line,*) 'wssmp_call: Using ',nthreads,' threads.'
            call C_OUT(2,0,1,line)
         endif
         call WSETMAXTHRDS(nthreads)
#ifdef SIZEOF_INT_P_IS_8
C
C     Reduced stack size to limit memory consumption
C
         if (nthreads .gt. 16) then
            call wsetmaxstack(0.25d0)
         else if (nthreads .gt. 8) then
            call wsetmaxstack(0.50d0)
         end if
#endif
C
C     Get default values for WSSMP
C
         IPARM(1) = 0
         IPARM(2) = 0
         IPARM(3) = 0
         call WSSMP( NTOT, IS(P_IA+1), IS(P_JA+1), dummy, dummy,
     1        IS(P_PERM+1), IS(P_INVP+1), dummy, 0, 0, dummy, 0, idummy,
     2        IPARM, DPARM)
C
C     Set them to what we need lateron
C
         IPARM(15) = 0 ! no restrictions on pivoting
         IPARM(16) = 1 ! slowest but best ordering
         IPARM(17) = 0 !
         IPARM(18) = 1 !
         IPARM(19) = 0 !
         IPARM(20) = 1 ! for ordering in IP methods?
         IPARM(31) = 1 ! want L D L^T decomposition
         IPARM(32) = 1 ! need D to see where first negative eigenvalue occurs

         IPARM(11) = 1
         DPARM(11) = 4.d-16
         DPARM(21) = 1.d+200

       endif
C
C     Now call WSSMP to perform ordering and symbolic factorization
C
      IPARM(2) = 1
      IPARM(3) = 2

      if( dump.ne.0 ) then
         open(45,file='matrix_unscaled.dat')
         write(45,*) 'Number of variables:'
         write(45,9901) NTOT
         write(45,*) 'IA'
         do i = 1, NTOT+1
            write(45,9901) IS(P_IA+i)
         enddo
         write(45,*) 'JA'
         do i = 1, IS(P_IA+NTOT+1)-1
            write(45,9901) IS(P_JA+i)
         enddo
      endif
      if( dump.eq.1 ) then
         write(45,*) 'IPARMS'
         do i = 1, 64
            write(45,9901) IPARM(i),i
         enddo
         write(45,*) 'DPARMS'
         do i = 1, 64
            write(45,9910) DPARM(i),i
         enddo
      endif
      call WSSMP(NTOT, IS(P_IA+1), IS(P_JA+1), dummy, dummy,
     1     IS(P_PERM+1), IS(P_INVP+1), dummy, 0, 0, dummy, 0, idummy,
     1     IPARM, DPARM)
      if( IPARM(64).ne.0 ) then
         write(line,*)
     1        'wssmp_call:  Error accurred in WSSMP-1-2 IPARM(64)=',
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
         write(line,*) 'wssmp_call:  Size of work space = ',IPARM(23)
         call C_OUT(1,0,1,line)
         write(line,*) 'wssmp_call:  Number of nonzeros in factor = ',
     1        IPARM(24)
         call C_OUT(1,0,1,line)
      endif
      if( dump.eq.1 ) then
         write(*,*) 'dump = ',dump
         stop
      endif

      LA = NZKKT
      LIWMA = 0

      goto 9999

 2000 continue
C
C     Do the factorization
C
      count = count + 1
      IPARM(2) = 3
      IPARM(3) = 3

      p_diag  = p_rwend
      p_avals = p_diag  + NTOT
      p_rwend = p_avals + NZCSC
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Copy elements into CSC format
C
      call DCOPY(NZCSC, 0.D0, 0, RW(p_avals+1), 1)
      do i = 1, NZKKT
         j = IS(P_PHSL2CSC+i)
         RW(p_avals+j) = RW(p_avals+j) + A(i)
      enddo
      if( .false. .and. count.eq.wantcount ) then
         do i = 1, NTOT
            do j = IS(P_IA+i), IS(P_IA+i+1)-1
               write(44,*) 'A(',i,',',IS(P_JA+j),')=',RW(p_avals+j)
            enddo
         enddo
      endif
C      if( dump.eq.2 ) then
C         do i = 1, 64
C            write(45,9901) IPARM(i)
C         enddo
C         do i = 1, 64
C            write(45,9910) DPARM(i)
C         enddo
C      endif
      if( (dump.eq.3 .and. count.eq.wantcount) .or. dump.eq.2 ) then
         write(45,*) 'AVALS for count ',count
         do i = 1, IS(P_IA+NTOT+1)-1
            write(45,9910) RW(p_avals+i)
         enddo
      endif
      if( dump.eq.2 ) then
         close(45)
         write(*,*) 'dump = ',dump
         IERR = 2155
         goto 9999
      endif
      call WSSMP(NTOT, IS(P_IA+1), IS(P_JA+1), RW(p_avals+1),
     1     RW(p_diag+1), IS(P_PERM+1), IS(P_INVP+1), dummy, 0, 0,
     1     dummy, 0, idummy, IPARM, DPARM)
      if( IPARM(64).ne.0 ) then
         write(line,*)
     1        'wssmp_call:  Error accurred in WSSMP-3 IPARM(64)=',
     2        IPARM(64)
         call C_OUT(2,0,1,line)
         if( IPARM(64).eq.-102 ) then
            IERR = 96
         else
            IERR = 592
         endif
         goto 9900
      endif
C
C     Determine number of negative eigenvalues
C
      NZERO = IPARM(21)
      NEIGS = 0
      do i = 1, NTOT
         if( RW(p_diag+i).lt.0.d0 ) then
            NEIGS = NEIGS + 1
         endif
      enddo
      L_INC_LA    = .false.
      L_INC_LIWMA = .false.

      goto 9999
C
C     Solve the system
C
 3000 continue
      IPARM(2) = 4
      IPARM(3) = 5
C      IPARM(3) = 4
      IPARM(6) = 0
C      IPARM(7) = 3
c      DPARM(6) = 1d-6
      if( dump.eq.3 .and. count.eq.wantcount ) then
         do j = 1, NEIGS
         write(45,*) 'RHS number ',j
         do i = 1, NTOT
            write(45,9910) RHS(i,j)
         enddo
         enddo
         close(45)
      endif
      call WSSMP(NTOT, IS(P_IA+1), IS(P_JA+1), dummy,
     1     dummy, IS(P_PERM+1), IS(P_INVP+1), RHS, NTOT, NEIGS,
     1     dummy, 0, idummy, IPARM, DPARM)
      if( IPARM(64).ne.0 ) then
         write(line,*)
     1        'wssmp_call:  Error accurred in WSSMP-4-5 IPARM(64)=',
     2        IPARM(64)
         call C_OUT(2,0,1,line)
         if( IPARM(64).eq.-102 ) then
            IERR = 96
         else
            IERR = 592
         endif
         goto 9900
      endif
      if( IPARM(6).gt.0 .and. QPRINT.ge.2 ) then
         write(line,*) 'Number of iterative refinement steps: ',IPARM(6)
         call C_OUT(1,0,1,line)
      endif
      if( dump.eq.3 .and. count.eq.wantcount ) then
         open(45,file='solution_unscaled.dat')
         do j = 1, NEIGS
         write(45,*) 'SOL number ',j
         do i = 1, NTOT
            write(45,9910) RHS(i,j)
         enddo
         enddo
         close(45)
         stop
      endif

      goto 9999
C
 9900 continue
      if( QPRINT.ge.1 ) then
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

      subroutine WSSMP_CALL_WS(TASK, N, M, NLB, NUB, NZA, NZKKT,
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
         LRW = 0
      elseif( TASK.eq.2 ) then
         LIW = 0
         LRW = N+M + NZKKT
      elseif( TASK.eq.3 ) then
         LIW = 0
         LRW = 0
      endif

      return
      end
