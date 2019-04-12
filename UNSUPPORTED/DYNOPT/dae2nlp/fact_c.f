C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine FACT_C(TASK, NZ, NY, NCOL, NAC, NELE, NNZFZZ, NNZFW,
     1     NNZAW, C, IRC, JCC, IELE_CRIT, LA48, FAST, IFACT, RFACT,
     2     LRW, RW, LIW, IW, IERR)
C
C     $Id: fact_c.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
C     Factorize overall basis matrix
C
C     TASK    1: Initialize data structures (based on element IELE_CRIT)
C                and factorize all elements
C             2: Factorize all elements
C     FAST    .false.: store all info for all elements
C             .true. : store IRN and KEEP only once (based on element
C                      IELE_CRIT)
C
      implicit none

      integer TASK, NZ, NY, NCOL, NAC, NELE, NNZFZZ, NNZFW, NNZAW
      double precision C(NELE*(NCOL*NNZFZZ+NNZFW+NNZAW))
      integer IRC(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer JCC(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer IELE_CRIT
      integer LA48
      logical FAST
      integer IFACT(*)          ! order: JCN(NNZ), IRN(LA48), KEEP(*)
                                ! (last 2 repeated if FAST)
      double precision RFACT(*)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR

      integer nkeep, lc, lfw, ljcn, lirn, lkeep, job, i
      integer lrfa, iele
      double precision dummy

      IERR  = 0
      nkeep = 7 + 10*(NCOL*(NZ+NY)+NAC)

      goto (1000, 2000) TASK

c      write(*,*) 'wrong TASK = ',TASK
      STOP
C
C     Initialize data structures of MA48
C
 1000 continue
      lc    = 1 + (IELE_CRIT-1)*(NCOL*NNZFZZ+NNZFW+NNZAW) + NCOL*NNZFZZ
      lfw   = 1 + NCOL*NNZFZZ
      ljcn  = 1
      lirn  = ljcn + NNZFW + NNZAW
      lkeep = lirn + LA48
      call FACT_FW(1, NZ, NY, NCOL, NAC, LA48, NNZFW+NNZAW, C(lc),
     1     IRC(lfw), JCC(lfw), dummy, IFACT(lirn), IFACT(ljcn),
     2     IFACT(lkeep), LRW, RW, LIW, IW, IERR)
      if( IERR.ne.0 ) goto 9999
C
C     Copy that information for every element
C
      if( .not.FAST ) then
         do i = 0, (NELE-1)*(LA48 + nkeep)-1
            IFACT(lirn + LA48 + nkeep + i) = IFACT(lirn+i)
         enddo
      endif

C
 2000 continue
C
C     First, factorize element IELE_CRIT
C
      lc    = 1 + (IELE_CRIT-1)*(NCOL*NNZFZZ+NNZFW+NNZAW) + NCOL*NNZFZZ
      lfw   = 1 + NCOL*NNZFZZ
      ljcn  = 1
      if( FAST ) then
         lirn = ljcn + NNZFW + NNZAW
      else
         lirn = ljcn + NNZFW + NNZAW + (IELE_CRIT-1)*(LA48+nkeep)
      endif
      lkeep = lirn + LA48
      lrfa  = 1    + (IELE_CRIT-1)*LA48
      call FACT_FW(2, NZ, NY, NCOL, NAC, LA48, NNZFW+NNZAW, C(lc),
     1     IRC(lfw), JCC(lfw), RFACT(lrfa), IFACT(lirn), IFACT(ljcn),
     2     IFACT(lkeep), LRW, RW, LIW, IW, IERR)
      if( IERR.ne.0 ) goto 9999
C
C     And now the others
C
      if( FAST ) then
         job = 3
      else
         job = 2
      endif
C
      do iele = 1, NELE
         if( iele.ne.IELE_CRIT ) then
            lc = 1 + (iele-1)*(NCOL*NNZFZZ+NNZFW+NNZAW) + NCOL*NNZFZZ
            if( .not.FAST ) then
               lirn  = ljcn + NNZFW + NNZAW + (iele-1)*(LA48+nkeep)
               lkeep = lirn + LA48
            endif
            lrfa = 1 + (iele-1)*LA48
            call FACT_FW(job, NZ, NY, NCOL, NAC, LA48, NNZFW+NNZAW,
     1           C(lc), IRC(lfw), JCC(lfw), RFACT(lrfa), IFACT(lirn),
     2           IFACT(ljcn), IFACT(lkeep), LRW, RW, LIW, IW, IERR)
            if( IERR.ne.0 ) goto 9999
         endif
      enddo

 9999 continue
      return
      end
