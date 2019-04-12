C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C ============================================================================
C
      SUBROUTINE APPSLN(T, NZ, NY, NU, NELE, NCOL, TI, X,
     1     ZVAL, DZVAL, YVAL, UVAL)
C !DEC$ ATTRIBUTES DLLEXPORT :: APPSLN
C
C     $Id: appsln.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Subroutine for evaluating the DAE variables at given time T
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
C     Version 1.0    ( 10-30-00 )
C
C     Description:
C     ============
C
C       evaluate Z, Y and U on time T based interpolation;
C          results in ZVAL, YVAL, UVAL
C
C     Parameters:
C
C     T        Input    time on which variable are to evaluated
C     NZ       Input    number of differential variables
C     NY       Input    number of algebraic variables
C     NU       Input    number of control variables
C     NELE     Input    number of elements
C     NCOL     Input    number of collocation points
C     TI       Input    boundaries of finite elements
C     X        Input    coefficients of collocation polynomials
C                         (for ordering see dae2nlp.f)
C     ZVAL     Output   value of differential variables at T
C     DZVAL    Output   value of derivative of differential variables at T
C     YVAL     Output   value of algebraic variables at T
C     UVAL     Output   value of control variables at T
C
      implicit none
      double precision T
      integer NZ, NY, NU, NELE, NCOL
      double precision TI(NELE+1)
      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ))
      double precision ZVAL(NZ)
      double precision DZVAL(NZ)
      double precision YVAL(NY)
      double precision UVAL(NU)

      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/

      integer iele, il, iu, lw, lu
      double precision theta(7), omega(7), dummy
C
C     Find element in which T is located (binary search)
C
      if( TI(2).ge.T ) then
         iele = 1
      elseif( TI(NELE).lt.T ) then
         iele = NELE
      else
         il = 2
         iu = NELE
 100     continue
         if( iu-il.gt.1 ) then
            iele = il+(iu-il)/2
            if( TI(iele).lt.T ) then
               il = IELE
            else
               iu = iele
            endif
            goto 100
         endif
         iele = il
      endif
C
C     Compute values
C
      lw = 1 + (iele-1)*(NCOL*(NZ+NY)+NZ)
      lu = 1 + NELE*(NCOL*(NZ+NY)+NZ) + NZ + (iele-1)*NCOL*NU
      call APPROX (3, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu), TI,
     1     iele, T, omega, theta, ZVAL, YVAL, UVAL)
      call APPROX (4, NZ, NY, NU, NELE, NCOL, COEF, X(lw), X(lu), TI,
     1     iele, T, omega, theta, DZVAL, dummy ,dummy)
C
C     That's it
C
      return
      end
