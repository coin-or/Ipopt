C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MULT_C(TRANS, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFZZ, NNZFW, NNZAW, NNZCW, NNZCWW, IRC, JCC, C, VIN, VOUT)
C
C     $Id: mult_c.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Multiply vector with overall basis matrix
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
C      integer IRC(NCOL*NNZFZZ+NNZFW+NNZAW+NNZCW)
C      integer JCC(NCOL*NNZFZZ+NNZFW+NNZAW+NNZCW)
      integer IRC(*)
      integer JCC(*)
      double precision C(NELE*(NCOL*NNZFZZ+NNZFW+NNZAW))
      double precision VIN(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
      double precision VOUT(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))

      integer lfz, lfw, lcw, lout, lc, iele
      double precision h

      lfz = 1
      lfw = lfz + NCOL*NNZFZZ
      lcw = lfw + NNZFW + NNZAW

      if( TRANS(1:1).eq.'t' .or. TRANS(1:1).eq.'T' ) then

         lout = 1
         lc   = 1
         do iele = 1, NELE
C
C     row with FZ^T
C
            call DCOPY(NZ, VIN(lout), 1, VOUT(lout), 1)
            call MULTMAT(TRANS, .true., NCOL*NNZFZZ, IRC(lfz), JCC(lfz),
     1           C(lc), VIN(lout+NZ), VOUT(lout))
            call DAXPY(NZ, -1d0, VIN(lout+NZ+NCOL*(NZ+NY)+NAC), 1,
     1           VOUT(lout), 1)
            lc = lc + NCOL*NNZFZZ
C
C     FW^T  and  AW^T
C
            call DCOPY(NCOL*(NZ+NY)+NAC, 0d0, 0, VOUT(lout+NZ), 1)
            call MULTMAT(TRANS, .true., NNZFW+NNZAW, IRC(lfw), JCC(lfw),
     1           C(lc), VIN(lout+NZ), VOUT(lout+NZ))
            lc = lc + NNZFW + NNZAW
C
C     CW^T
C
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .true., NZ, NY, NU, NP, NCOL, h,
     1           NNZCW, NNZCWW, IRC(lcw), JCC(lcw),
     1           VIN(lout+NZ+NCOL*(NZ+NY)+NAC), VOUT(lout+NZ))
C
            lout = lout + NZ + NCOL*(NZ+NY) + NAC

         enddo

         call DCOPY(NZ, VIN(lout), 1, VOUT(lout), 1)

      else
C
C     start with the initial conditions
C
         call DCOPY(NZ, VIN, 1, VOUT, 1)
C
C     loop over all elements
C
         lout = 1
         lc   = 1
         do iele = 1, NELE
C
C     FZ
C
            call DCOPY(NCOL*(NZ+NY)+NAC, 0d0, 0, VOUT(lout+NZ), 1)
            call MULTMAT(TRANS, .true., NCOL*NNZFZZ, IRC(lfz), JCC(lfz),
     1           C(lc), VIN(lout), VOUT(lout+NZ))
            lc = lc + NCOL*NNZFZZ
C
C     FW  and  AW
C
            call MULTMAT(TRANS, .true., NNZFW+NNZAW, IRC(lfw), JCC(lfw),
     1           C(lc), VIN(lout+NZ), VOUT(lout+NZ))
            lc = lc + NNZFW + NNZAW
C
C     Continuity equations
C
            call DCOPY(NZ, VIN(lout+NZ+NCOL*(NZ+NY)+NAC), 1,
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC), 1)
            call DAXPY(NZ, -1d0, VIN(lout), 1,
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC), 1)
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .true., NZ, NY, NU, NP, NCOL, h,
     1           NNZCW, NNZCWW, IRC(lcw), JCC(lcw), VIN(lout+NZ),
     1           VOUT(lout+NZ+NCOL*(NZ+NY)+NAC))
C
            lout = lout + NZ + NCOL*(NZ+NY) + NAC

         enddo

      endif

      return
      end
