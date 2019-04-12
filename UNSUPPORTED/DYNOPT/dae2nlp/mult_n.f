C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine MULT_N(TRANS, NZ, NY, NU, NP, NCOL, NAC, NELE, TI,
     1     NNZFU, NNZAU, NNZFPP, NNZCU, NNZCUU, IRN, JCN, N, VIN, VOUT)
C
C     $Id: mult_n.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Multiply vector with overall nonbasis matrix
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
      integer NNZFU
      integer NNZAU
      integer NNZFPP
      integer NNZCU
      integer NNZCUU(NCOL)
C      integer IRN(NNZFU+NNZAU+NCOL*NNZFPP+NNZCU)
C      integer JCN(NNZFU+NNZAU+NCOL*NNZFPP+NNZCU)
      integer IRN(*)
      integer JCN(*)
C      double precision N(NELE*(NNZFU+NNZAU+NCOL*NNZFPP))
      double precision N(*)
      double precision VIN(*)
      double precision VOUT(*)

      integer lfu, lfp, lcu, lout, lin, ln, lu, lp, iele
      double precision h

      lfu = 1
      lfp = lfu + NNZFU + NNZAU
      lcu = lfp + NCOL*NNZFPP

      if( TRANS(1:1).eq.'t' .or. TRANS(1:1).eq.'T' ) then

         lin = 1 + NZ
         lu  = 1
         lp  = 1 + NELE*(NCOL*NU-NAC) 
         ln  = 1

         call DCOPY(NELE*(NCOL*NU-NAC) + NP, 0d0, 0, VOUT, 1)

         do iele = 1, NELE
C
C     FU^T  and  AU^T
C
            call MULTMAT(TRANS, .true., NNZFU+NNZAU, IRN(lfu), JCN(lfu),
     1           N(ln), VIN(lin), VOUT(lu))
            ln  = ln  + NNZFU + NNZAU
C
C     FP^T
C
            call MULTMAT(TRANS, .true., NCOL*NNZFPP, IRN(lfp), JCN(lfp),
     1           N(ln), VIN(lin), VOUT(lp))
            ln  = ln  + NCOL*NNZFPP
            lin = lin + NCOL*(NZ+NY) + NAC
C
C     CU^T
C
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .true., NZ, NY, NU, NP, NCOL, h,
     1           NNZCU, NNZCUU, IRN(lcu), JCN(lcu), VIN(lin), VOUT(lu))
C
            lin = lin + NZ
            lu  = lu  + NCOL*NU - NAC

         enddo

      else

         call DCOPY(NZ + NELE*(NCOL*(NZ+NY) + NZ + NAC), 0d0, 0,
     1        VOUT, 1)
C
C     loop over all elements
C
         lout = 1 + NZ
         lu   = 1
         lp   = 1 + NELE*(NCOL*NU-NAC)
         ln   = 1

         do iele = 1, NELE
C
C     FU  and  AU
C
            call MULTMAT(TRANS, .true., NNZFU+NNZAU, IRN(lfu), JCN(lfu),
     1           N(ln), VIN(lu), VOUT(lout))
            ln   = ln   + NNZFU + NNZAU
C
C     FP
C
            call MULTMAT(TRANS, .true., NCOL*NNZFPP, IRN(lfp), JCN(lfp),
     1           N(ln), VIN(lp), VOUT(lout))
            ln   = ln   + NCOL*NNZFPP
            lout = lout + NCOL*(NZ+NY) + NAC
C
C     Continuity equations
C
            h = TI(iele+1) - TI(iele)
            call MULTCW(TRANS, .true., NZ, NY, NU, NP, NCOL, h,
     1           NNZCU, NNZCUU, IRN(lcu), JCN(lcu), VIN(lu), VOUT(lout))
C
            lout = lout + NZ
            lu   = lu   + NCOL*NU - NAC

         enddo

      endif

      return
      end
