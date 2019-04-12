C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine SOLVE_C(TRANS, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFZZ, NNZFW, NNZAW, NNZCW, NNZCWW, IRC, JCC, C, LA48, FAST,
     2     IFACT, RFACT, VIN, VOUT, rhs, LRW, RW, LIW, IW, IERR)
C
C     $Id: solve_c.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Solve system involving overall basis matrix
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none

      character*(*) TRANS       !if = 't' or 'T', compute vout += CW^T*vin
                                !otherwise vout += CW*vin
                                !NOTE: VOUT is NOT INITIALIZED here!!!
      integer NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision TI(NELE+1)
      integer NNZFZZ
      integer NNZFW
      integer NNZAW
      integer NNZCW
      integer NNZCWW(NCOL)
      integer IRC(NCOL*NNZFZZ+NNZFW+NNZAW+NNZCW)
      integer JCC(NCOL*NNZFZZ+NNZFW+NNZAW+NNZCW)
      double precision C(NELE*(NCOL*NNZFZZ+NNZFW+NNZAW))
      integer LA48
      logical FAST
      integer IFACT(*)
      double precision RFACT(*)
      double precision VIN(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
      double precision VOUT(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
      double precision rhs(NCOL*(NZ+NY)+NAC) ! work space
      integer LRW
      double precision RW(LRW)  !real workspace for SOLVE_FW
      integer LIW
      integer IW(LIW)           !integer workspace for SOLVE_FW
      integer IERR              !error from SOLVE_FW

      integer lfz, lfw, lcw, lout, lc, iele, nkeep, lrfa, lirn
      integer lkeep
      double precision h

      IERR  = 0
      nkeep = 7 + 10*(NCOL*(NZ+NY)+NAC)
      lfz   = 1
      lfw   = lfz + NCOL*NNZFZZ
      lcw   = lfw + NNZFW + NNZAW

      if( TRANS(1:1).eq.'t' .or. TRANS(1:1).eq.'T' ) then

         lout = 1 + (NELE-1)*(NZ + NCOL*(NZ+NY)+NAC)
         lc   = 1 + (NELE-1)*(NCOL*NNZFZZ + NNZFW + NNZAW)
         lrfa = 1 + (NELE-1)*LA48
         if( FAST ) then
            lirn = 1 + NNZFW + NNZAW
         else
            lirn = 1 + NNZFW + NNZAW + (NELE-1)*(LA48 + nkeep)
         endif
         lkeep = lirn + LA48
C
C     start with very last identity block
C
         call DCOPY(NZ, VIN(lout + NZ + NCOL*(NZ+NY) + NAC), 1,
     1        VOUT(lout + NZ + NCOL*(NZ+NY) + NAC), 1)

         do iele = NELE, 1, -1
C
C     row with  FW^T, AW^T  and  CW^T
C
            call DCOPY(NCOL*(NZ+NY)+NAC, VIN(lout+NZ), 1, rhs, 1)
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .false., NZ, NY, NU, NP, NCOL, h,
     1           NNZCW, NNZCWW, IRC(lcw), JCC(lcw),
     1           VOUT(lout + NZ + NCOL*(NZ+NY) + NAC), rhs)
            call SOLVE_FW(TRANS, NZ, NY, NCOL, NAC, LA48, NNZFW+NNZAW,
     1           RFACT(lrfa), IFACT(lirn), IFACT(lkeep), rhs,
     2           VOUT(lout+NZ), LRW, RW, LIW, IW, IERR)
            if( IERR.ne.0 ) goto 9999
C
C     row with FZ^T
C
            call DCOPY(NZ, VIN(lout), 1, VOUT(lout), 1)
            call DAXPY(NZ, 1d0, VOUT(lout + NZ + NCOL*(NZ+NY) + NAC), 1,
     1           VOUT(lout), 1)
            call MULTMAT(TRANS, .false., NCOL*NNZFZZ, IRC(lfz),
     1           JCC(lfz), C(lc), VOUT(lout+NZ), VOUT(lout))
C
            lout = lout - (NZ + NCOL*(NZ+NY) + NAC)
            lc   = lc   - (NCOL*NNZFZZ + NNZFW + NNZAW)
            lrfa = lrfa - LA48
            if( .not.FAST ) then
               lirn  = lirn  - (LA48 + nkeep)
               lkeep = lkeep - (LA48 + nkeep)
            endif

         enddo

      else
C
C     start with the initial conditions
C
         call DCOPY(NZ, VIN, 1, VOUT, 1)
C
C     loop over all elements
C
         lout  = 1
         lc    = 1
         lrfa  = 1
         lirn  = 1 + NNZFW + NNZAW
         lkeep = lirn + LA48
         do iele = 1, NELE
C
C     row with FZ and FW
C
            call DCOPY(NCOL*(NZ+NY)+NAC, VIN(lout+NZ), 1, rhs, 1)
C
C     FZ
C
            call MULTMAT(TRANS, .false., NCOL*NNZFZZ, IRC(lfz),
     1           JCC(lfz), C(lc), VOUT(lout), rhs)
C
C     FW  and  AW
C
            call SOLVE_FW(TRANS, NZ, NY, NCOL, NAC, LA48, NNZFW+NNZAW,
     1           RFACT(lrfa), IFACT(lirn), IFACT(lkeep), rhs,
     2           VOUT(lout+NZ), LRW, RW, LIW, IW, IERR)
            if( IERR.ne.0 ) goto 9999
C
C     row with CW in it
C
            call DCOPY(NZ, VIN(lout+NZ+NCOL*(NZ+NY)+NAC), 1,
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC), 1)
            call DAXPY(NZ, 1d0, VOUT(lout), 1,
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC), 1)
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .false., NZ, NY, NU, NP, NCOL, h,
     1           NNZCW, NNZCWW, IRC(lcw), JCC(lcw), VOUT(lout+NZ),
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC))
C
            lout = lout + NZ + NCOL*(NZ+NY) + NAC
            lc   = lc   + NCOL*NNZFZZ + NNZFW + NNZAW
            lrfa = lrfa + LA48
            if( .not.FAST ) then
               lirn  = lirn  + LA48 + nkeep
               lkeep = lkeep + LA48 + nkeep
            endif

         enddo

      endif

 9999 continue
      return
      end
