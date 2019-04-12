C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine FILLJAC(NZ, NY, NU, NP, NCOL, NAC, NELE, TI, X, IVARC,
     1     IVARC1, NNZFZZ, NNZFW, NNZFWZ, NNZFWD, NNZFWY, NNZFWU, NNZFU,
     2     NNZFUZ, NNZFUD, NNZFUY, NNZFUU, NNZFPP, NNZAW, NNZAU,
     3     IRC, JCC, IRN, JCN, C, N, zval, lddf, df,
     4     NNZAC, irac, jcac, aac)
C
C     $Id: filljac.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute nonzeros for overall Jacobian
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
C      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision X(*)
      integer IVARC(NCOL*(NZ+NY+NU))
      integer IVARC1(NCOL*(NZ+NY+NU)) ! inverse of IVARC (only needed it NAC>0)
      integer NNZFZZ            ! nnz in FZ per col. eqn.
      integer NNZFW             ! total nnz in FW
      integer NNZFWZ(NCOL)      ! nnz in FW for dz per col. eqn.
      integer NNZFWD(NCOL)      ! nnz in FW for dz per col. eqn.
      integer NNZFWY(NCOL)      ! nnz in FW for y per col. eqn.
      integer NNZFWU(NCOL)      ! nnz in FW for u per col. eqn.
      integer NNZFU             ! total nnz in FU
      integer NNZFUZ(NCOL)      ! nnz in FU for dz per col. eqn.
      integer NNZFUD(NCOL)      ! nnz in FU for dz per col. eqn.
      integer NNZFUY(NCOL)      ! nnz in FU for y per col. eqn.
      integer NNZFUU(NCOL)      ! nnz in FU for u per col. eqn.
      integer NNZFPP            ! nnz in FP per col. eqn.
      integer NNZAW
      integer NNZAU
      integer IRC(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer JCC(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer IRN(NNZFU+NCOL*NNZFPP+NNZAU)
      integer JCN(NNZFU+NCOL*NNZFPP+NNZAU)
      double precision C(NELE*(NCOL*NNZFZZ+NNZFW+NNZAW))
      double precision N(*)  !(NELE*(NNZFU+NCOL*NNZFPP+NNZAU))
      double precision zval(NZ) ! work space
      integer lddf              ! (at least NZ+NY)
      double precision df(lddf,2*NZ+NY+NU+NP)
      integer NNZAC
      integer irac(NNZAC)
      integer jcac(NNZAC)
      double precision aac(NNZAC)

      integer iele
      integer lc, ln, lw, lu, lp, i, l

      lc = 1
      ln = 1
      lw = 1
      lu = NZ + NELE*(NCOL*(NZ+NY)+NZ) + 1
      lp = lu + NELE*(NCOL*NU)
      do iele = 1, NELE
C
C     Evaluate the Jacobian for that element
C
         call FILLBLOK(NZ, NY, NU, NP, NCOL, NAC, NELE, X(lw), X(lu),
     1        X(lp), TI, iele, IVARC, IVARC1, NNZFZZ, NNZFW, NNZFWZ,
     2        NNZFWD, NNZFWY, NNZFWU, NNZFU, NNZFUZ, NNZFUD, NNZFUY,
     3        NNZFUU, NNZFPP, NNZAW, NNZAU, IRC, JCC,
     4        C(lc), IRN, JCN, N(ln), lddf, df, zval,
     5        NNZAC, irac, jcac, aac)
C
C     Increase pointers
C
         lc = lc + NCOL*NNZFZZ + NNZFW + NNZAW
         ln = ln + NNZFU + NCOL*NNZFPP + NNZAU
         lw = lw + NCOL*(NZ+NY) + NZ
         lu = lu + NCOL*NU
C
      enddo

 9999 continue
      return
      end
