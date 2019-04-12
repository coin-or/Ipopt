C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine COLRHS(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, ZINIT, X,
     1     RHS, zval)
C
C     $Id: colrhs.f 531 2004-03-11 01:31:07Z andreasw $
C
C     compute right hand side for the collocation system
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none
      integer NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision TI(NELE+1)
      double precision ZINIT(NZ)
C      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision X(*)
      double precision RHS(NZ + NELE*(NCOL*(NZ+NY)+NZ+NAC))
      double precision zval(NZ) ! work space

      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/

      integer lw, lu, lp, lrhs, iele, i, l, k
      double precision h, hrho, trho, dummy

      lw   = 1
      lu   = NZ + NELE*(NCOL*(NZ+NY)+NZ) + 1
      lp   = lu + NELE*NCOL*NU
      lrhs = 1
C
C     Initial conditions
C
      call DCOPY(NZ, X, 1, RHS(lrhs), 1)
      call DAXPY(NZ, -1d0, ZINIT, 1, RHS(lrhs), 1)
      lrhs = lrhs + NZ
C
C     Loop over all elements
C
      do iele = 1, NELE
C
C     Now loop over the collocation equations
C
         h = TI(iele+1) - TI(iele)
         do k = 1, NCOL
C
C     Compute values of Z at that collocation point
C
            hrho = h * RHO(k)
            trho = TI(iele) + hrho
            call APPROX(1, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu),
     1           TI, iele, trho, OMEGAQ(1,k), dummy, zval, dummy, dummy)
C
C     Compute DAE equations from model
C
            call DAEMODEL_F(NZ, NY, NU, NP, trho,
     1           zval, X(lw+NZ+(k-1)*(NZ+NY)),
     1           X(lw+NZ+(k-1)*(NZ+NY)+NZ), X(lu+(k-1)*NU), X(lp),
     2           RHS(lrhs))
C
            lrhs = lrhs + NZ+NY

         enddo
C
C     Additional constraints
C
         if( NAC.gt.0 ) then
            call ADDCON_F(NZ, NY, NU, NCOL, NAC, X(lw+NZ), X(lu),
     1           RHS(lrhs))
            lrhs = lrhs + NAC
         endif
C
C     Continuity equations
C
         trho = TI(iele+1)
         call APPROX(1, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu), TI,
     1           iele, trho, OMEGA1, dummy, zval, dummy, dummy)
         call DCOPY(NZ, X(lw+NZ+NCOL*(NZ+NY)), 1, RHS(lrhs), 1)
         call DAXPY(NZ, -1d0, zval, 1, RHS(lrhs), 1)
         lrhs = lrhs + NZ
C
         lw = lw + NCOL*(NZ+NY) + NZ
         lu = lu + NCOL*NU

      enddo

      return
      end
