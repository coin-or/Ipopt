C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C
C     $Id: addcon.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      subroutine ADDCON_INIT(NZ, NY, NU, NCOL, NDEGU, NAC)
C !DEC$ ATTRIBUTES DLLEXPORT :: ADDCON_INIT
      implicit none
      integer NZ, NY, NU, NCOL, NDEGU, NAC
      integer ADDCONFLAG
      common /ADDCON/ ADDCONFLAG
      save /ADDCON/
      if( NDEGU.lt.0 .or. NDEGU.ge.NCOL-1 ) then
         ADDCONFLAG = 0
         NAC = 0
      else
         ADDCONFLAG = 10*NCOL+NDEGU
         NAC = NU*(NCOL-NDEGU-1)
      endif
      return
      end

      subroutine ADDCON_NNZ(NZ, NY, NU, NNZAC)
      implicit none
      integer NZ, NY, NU, NNZAC
      integer ADDCONFLAG
      common /ADDCON/ ADDCONFLAG
      save /ADDCON/
      if( ADDCONFLAG.eq.0 ) then
         NNZAC = 0
      elseif( ADDCONFLAG.eq.20 ) then
         NNZAC = 2*NU
      elseif( ADDCONFLAG.eq.31 ) then
         NNZAC = 3*NU
      elseif( ADDCONFLAG.eq.30 ) then
         NNZAC = 4*NU
      elseif( ADDCONFLAG.eq.40 ) then
         NNZAC = 6*NU
      elseif( ADDCONFLAG.eq.50 ) then
         NNZAC = 8*NU
      else
         write(*,*) 'ADDCON_NNZ: Invalid ADDCONFLAG = ',ADDCONFLAG
         stop
      endif
      return
      end

      subroutine ADDCON_F(NZ, NY, NU, NCOL, NAC, W, U, ADDCON)
      implicit none
      integer NZ, NY, NU, NCOL, NAC
      double precision W(NCOL*(NZ+NY))
      double precision U(NCOL*NU)
      double precision ADDCON(NAC)
      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/
      integer ADDCONFLAG
      common /ADDCON/ ADDCONFLAG
      save /ADDCON/
      integer i
      if( ADDCONFLAG.eq.20 ) then
         do i = 1, NU
            ADDCON(i) = U(i) - U(NU+i)
         enddo
      elseif( ADDCONFLAG.eq.31 ) then
         do i = 1, NU
            ADDCON(i) = (RHO(2)-RHO(3))*U(i) + (RHO(3)-RHO(1))*U(NU+i) +
     1           (RHO(1)-RHO(2))*U(2*NU+i)
         enddo
      elseif( ADDCONFLAG.eq.30 ) then
         do i = 1, NU
            ADDCON(   i) = U(   i) - U(  NU+i)
            ADDCON(NU+i) = U(NU+i) - U(2*NU+i)
         enddo
      elseif( ADDCONFLAG.eq.40 ) then
         do i = 1, NU
            ADDCON(     i) = U(     i) - U(  NU+i)
            ADDCON(  NU+i) = U(  NU+i) - U(2*NU+i)
            ADDCON(2*NU+i) = U(2*NU+i) - U(3*NU+i)
         enddo
      elseif( ADDCONFLAG.eq.50 ) then
         do i = 1, NU
            ADDCON(     i) = U(     i) - U(  NU+i)
            ADDCON(  NU+i) = U(  NU+i) - U(2*NU+i)
            ADDCON(2*NU+i) = U(2*NU+i) - U(3*NU+i)
            ADDCON(3*NU+i) = U(3*NU+i) - U(4*NU+i)
         enddo
      elseif( ADDCONFLAG.ne.0 ) then
         write(*,*) 'ADDCON_F: Invalid ADDCONFLAG = ',ADDCONFLAG
         stop
      endif
      return
      end

      subroutine ADDCON_DF(NZ, NY, NU, NCOL, NAC, W, U, NNZAC, IRAC,
     1     JCAC, AAC)
      implicit none
      integer NZ, NY, NU, NCOL, NAC
      double precision W(NCOL*(NZ+NY))
      double precision U(NCOL*NU)
      integer NNZAC
      integer IRAC(NNZAC)
      integer JCAC(NNZAC)
      double precision AAC(NNZAC)
      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/
      integer ADDCONFLAG
      common /ADDCON/ ADDCONFLAG
      save /ADDCON/
      integer lu, i
      lu = NCOL*(NZ+NY)
      if( ADDCONFLAG.eq.20 ) then
         do i = 1, NU
            IRAC(   i) = i
            JCAC(   i) = lu + i
            AAC (   i) = 1.d0
            IRAC(NU+i) = i
            JCAC(NU+i) = lu + NU + i
            AAC (NU+i) = -1.d0
         enddo
      elseif( ADDCONFLAG.eq.31 ) then
         do i = 1, NU
            IRAC(     i) = i
            JCAC(     i) = lu + i
            AAC (     i) = RHO(2)-RHO(3)
            IRAC(  NU+i) = i
            JCAC(  NU+i) = lu + NU + i
            AAC (  NU+i) = RHO(3)-RHO(1)
            IRAC(2*NU+i) = i
            JCAC(2*NU+i) = lu + 2*NU + i
            AAC (2*NU+i) = RHO(1)-RHO(2)
         enddo
      elseif( ADDCONFLAG.eq.30 ) then
         do i = 1, NU
            IRAC(     i) = i
            JCAC(     i) = lu + i
            AAC (     i) = 1.d0
            IRAC(  NU+i) = i
            JCAC(  NU+i) = lu + NU + i
            AAC (  NU+i) = -1.d0
            IRAC(2*NU+i) = NU+i
            JCAC(2*NU+i) = lu + NU + i
            AAC (2*NU+i) = 1.d0
            IRAC(3*NU+i) = NU+i
            JCAC(3*NU+i) = lu + 2*NU + i
            AAC (3*NU+i) = -1.d0
         enddo
      elseif( ADDCONFLAG.eq.40 ) then
         do i = 1, NU
            IRAC(     i) = i
            JCAC(     i) = lu + i
            AAC (     i) = 1.d0
            IRAC(  NU+i) = i
            JCAC(  NU+i) = lu + NU + i
            AAC (  NU+i) = -1.d0
            IRAC(2*NU+i) = NU+i
            JCAC(2*NU+i) = lu + NU + i
            AAC (2*NU+i) = 1.d0
            IRAC(3*NU+i) = NU+i
            JCAC(3*NU+i) = lu + 2*NU + i
            AAC (3*NU+i) = -1.d0
            IRAC(4*NU+i) = 2*NU+i
            JCAC(4*NU+i) = lu + 2*NU + i
            AAC (4*NU+i) = 1.d0
            IRAC(5*NU+i) = 2*NU+i
            JCAC(5*NU+i) = lu + 3*NU + i
            AAC (5*NU+i) = -1.d0
         enddo
      elseif( ADDCONFLAG.eq.50 ) then
         do i = 1, NU
            IRAC(     i) = i
            JCAC(     i) = lu + i
            AAC (     i) = 1.d0
            IRAC(  NU+i) = i
            JCAC(  NU+i) = lu + NU + i
            AAC (  NU+i) = -1.d0
            IRAC(2*NU+i) = NU+i
            JCAC(2*NU+i) = lu + NU + i
            AAC (2*NU+i) = 1.d0
            IRAC(3*NU+i) = NU+i
            JCAC(3*NU+i) = lu + 2*NU + i
            AAC (3*NU+i) = -1.d0
            IRAC(4*NU+i) = 2*NU+i
            JCAC(4*NU+i) = lu + 2*NU + i
            AAC (4*NU+i) = 1.d0
            IRAC(5*NU+i) = 2*NU+i
            JCAC(5*NU+i) = lu + 3*NU + i
            AAC (5*NU+i) = -1.d0
            IRAC(6*NU+i) = 3*NU+i
            JCAC(6*NU+i) = lu + 3*NU + i
            AAC (6*NU+i) = 1.d0
            IRAC(7*NU+i) = 3*NU+i
            JCAC(7*NU+i) = lu + 4*NU + i
            AAC (7*NU+i) = -1.d0
         enddo
      elseif( ADDCONFLAG.ne.0 ) then
         write(*,*) 'ADDCON_F: Invalid ADDCONFLAG = ',ADDCONFLAG
         stop
      endif
      return
      end
