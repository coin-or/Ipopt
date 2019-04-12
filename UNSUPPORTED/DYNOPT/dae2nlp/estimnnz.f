C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine ESTIMNNZ(NZ, NY, NU, NP, NCOL, NAC, LDDF, DF,
     1     NNZA, NNZC, NNZAC, IRAC, JCAC, AAC, WDUMMY, UDUMMY, IW)
C
C     $Id: estimnnz.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute demand of space to store A, and upper bound of nonzeros in C
C     for any partition
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none

      integer NZ, NY, NU, NP, NCOL, NAC
      integer LDDF
      double precision DF(LDDF, 2*NZ+NY+NU+NP)
      integer NNZA, NNZC
      integer NNZAC
      integer IRAC(NNZAC)
      integer JCAC(NNZAC)
      double precision AAC(NNZAC)
      double precision WDUMMY(NCOL*(NZ+NY))
      double precision UDUMMY(NCOL*NU)
      integer IW(NZ + NCOL*(NZ+NY+NU) + NP)

      integer i, j, k, liw
C
C     Shouldn't be necessary, but initialize to zero
C
      do i = 1, 2*NZ+NY+NU+NP
         call DCOPY(NZ+NY, 0d0, 0, DF(1,i), 1)
      enddo
      call DAEMODEL_DF_STRUC(NZ, NY, NU, NP, DF, LDDF)
C
C     count nonzeros for Z's
C
      do j = 1, NZ
         IW(j) = 0
         do i = 1, NZ+NY
            if( DF(i,j).ne.0d0 ) IW(j) = IW(j) + NCOL
         enddo
      enddo
C
C     nonzeros for the Zdot's
C
      do j = 1, NZ
         IW(NZ+j) = IW(j)
         do i = 1, NZ+NY
            if( DF(i,NZ+j).ne.0d0 ) IW(NZ+j) = IW(NZ+j) + 1
         enddo
      enddo
C
C     nonzeros for y's
C
      liw = 2*NZ
      do j = 1, NY
         IW(liw+j) = 0
         do i = 1, NZ+NY
            if( DF(i,2*NZ+j).ne.0d0 ) IW(liw+j) = IW(liw+j) + 1
         enddo
      enddo
C
C     copy for remaining collocation points
C
      liw = NZ
      do i = 2, NCOL
         do j = 1, NZ+NY
            IW(liw+NZ+NY+j) = IW(liw+j)
         enddo
         liw = liw + NZ+NY
      enddo
C
C     nonzeros for u's
C
      liw = NZ + NCOL*(NZ+NY)
      do j = 1, NU
         IW(liw+j) = 0
         do i = 1, NZ+NY
            if( DF(i,2*NZ+NY+j).ne.0d0 ) IW(liw+j) = IW(liw+j) + 1
         enddo
      enddo
C
C     copy for remaining collocation points
C
      liw = NZ + NCOL*(NZ+NY)
      do i = 2, NCOL
         do j = 1, NU
            IW(liw+NU+j) = IW(liw+j)
         enddo
         liw = liw + NU
      enddo
C
C     finally the p's
C
      liw = NZ + NCOL*(NZ+NY+NU)
      do j = 1, NP
         IW(liw+j) = 0
         do i = 1, NZ+NY
            if( DF(i,2*NZ+NY+NU+j).ne.0d0 ) IW(liw+j) = IW(liw+j) + NCOL
         enddo
      enddo
C
C     Well, we are not done yet :)  Let's look at the additional constraints
C
      if( NAC.gt.0 ) then
         call ADDCON_DF(NZ, NY, NU, NCOL, NAC, WDUMMY, UDUMMY, NNZAC,
     1        IRAC, JCAC, AAC)
         do i = 1, NNZAC
            j = NZ+JCAC(i)
            IW(j) = IW(j) + 1
         enddo
      endif
C
C     OK, now we can count the nonzeros
C
      NNZA = 0
      do i = 1, NZ + NCOL*(NZ+NY+NU) + NP
         NNZA = NNZA + IW(i)
      enddo
C
C     Now sort the number of nonzeros in [FW FU] into descending order
C
      do i = 1, NCOL*(NZ+NY+NU)
         do j = i+1, NCOL*(NZ+NY+NU)
            if( IW(NZ+i).lt.IW(NZ+j) ) then
               k = IW(NZ+i)
               IW(NZ+i) = IW(NZ+j)
               IW(NZ+j) = k
            endif
         enddo
      enddo
C
C     Instead of using the inefficient version above, you may choose to use
C     another method like quick sort, e.g. implemented in Harwell's KB06AI
C
C      call KB06AI(IW(NZ+1), NCOL*(NZ+NY+NU))
C
C     The worst estimate for C is simply the sum of the first largest
C     column counts
C
      NNZC = 0
      do i = NZ + 1, NZ + NCOL*(NZ+NY) + NAC
         NNZC = NNZC + IW(i)
      enddo
C
C     The End
C
      return
      end
