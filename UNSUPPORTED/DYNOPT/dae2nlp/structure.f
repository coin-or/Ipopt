C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine STRUCTURE(NZ, NY, NU, NP, NCOL, NAC,
     1     IVARC1, NNZA, IR, JC, LDDF, DF, NNZAC, IRAC, JCAC, AAC,
     2     NNZFZZ, NNZFW, NNZFWZ, NNZFWD, NNZFWY, NNZFWU, NNZAW, NNZCW,
     3     NNZCWW, NNZFU, NNZFUZ, NNZFUD, NNZFUY, NNZFUU, NNZFPP, NNZAU,
     4     NNZCU, NNZCUU, WDUMMY, UDUMMY)
C
C     $Id: structure.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Initialize structure of elemental Jacobian
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
C     IVARC1  assume, xc = (dz1, y1, dz2, y2, ..., dzN, yN, u1, ..., uN)
C             then IVARC1(i) is column number in C corresponding to xc(i)
C             (if <= NCOL*(NZ+NY) ), and IVARC1(i)-NCOL*(NZ+NY) is col no.
C             in N
C
C     Structure of Jacobian for one element (after reordering)
C     (you will need a REALLY broad terminal to see that one... :)
C
C
C        |---------------------------------   C   ------------------------------------...|           |...-------------------------------  N  --------------------------------|
C

C      FZ(1)   FWZ(1,1)+FWD(1)   FWY(1)   FWU(1)   ...      FWZ(1,K)         0        0       FUZ(1,1)+FWD(1)   FUY(1)   FUU(1)  ...      FUZ(1,K)         0        0      FP(1)
C      FZ(2)       FWZ(2,1)        0        0      ...      FWZ(2,K)         0        0           FUZ(2,1)        0        0     ...      FUZ(2,K)         0        0      FP(2)
C        .            .            .        .                  .             .        .              .            .        .                 .             .        .        .
C        .            .            .        .                  .             .        .              .            .        .                 .             .        .        .
C        .            .            .        .                  .             .        .              .            .        .                 .             .        .        .
C      FZ(K)       FWZ(K,1)        0        0      ...   FWZ(K,K)+FWD(K)   FWY(K)   FWU(K)        FUZ(K,1)        0        0     ...   FUZ(K,K)+FUD(K)   FUY(K)   FUU(K)   FP(K)
C        <---------------------------------   AW  ------------------------------------>              <---------------------------   AU   ---------------------------->       0
C       -I          CW(1)          0        0      ...       CW(K)           0        0              CU(1)        0        0     ...       CU(K)           0        0        0       I
C
C
C     Order within JC and IR:  (offsets w.r.t. matrix above)
C
C     [  FZ(1)      FZ(2)      ...   FZ(K)
C            (column offset 0, row offset 0)
C        FWZ(1,1)   FWZ(1,2)   ...   FWZ(1,K)   FWD(1)   FWY(1)   FWU(1)
C          :
C          :
C        FWZ(K,1)   FWZ(K,2)   ...   FWZ(K,K)   FWD(K)   FWY(K)   FWU(K)
C            (column offset NZ, row offset 0)
C        AW
C            (column offset NZ, row offset 0)
C        CW(1)      CW(2)      ...   CW(K)
C            (column offset NZ, row offset NCOL*(NZ+NY)+NAC)
C
C        FUZ(1,1)   FUZ(1,2)   ...   FUZ(1,K)   FUD(1)   FUY(1)   FUU(1)
C          :
C          :
C        FUZ(K,1)   FUZ(K,2)   ...   FUZ(K,K)   FUD(K)   FUY(K)   FUU(K)
C            (column offset NZ+NCOL*(NZ+NY)+NAC, row offset 0)
C        AU
C            (column offset NZ+NCOL*(NZ+NY)+NAC, row offset 0)
C        FP(1)      FP(2)      ...   FP(K)
C            (column offset NZ+NCOL*(NZ+NY+NU)+NAC, row offset 0)
C        CU(1)      CU(2)      ...   CU(K)  ]
C            (column offset NZ+NCOL*(NZ+NY)+NAC, row offset NCOL*(NZ+NY)+NAC)
C
C      (Second part here has again column offset 0)
C
      implicit none
      integer NZ, NY, NU, NP, NCOL, NAC
      integer IVARC1(NCOL*(NZ+NY+NU))
      integer NNZA              ! number of nonzeros in A (from ESTIMNNZ)
                                ! (doesn't include entries for CW CU, but
                                ! for AW and AU)
      integer IR(NNZA+NCOL*NZ)
      integer JC(NNZA+NCOL*NZ)
      integer LDDF
      double precision DF(LDDF, 2*NZ+NY+NU+NP) ! work space
      integer NNZAC             ! total nnz for additional constraints
      integer IRAC(NNZAC)       ! work space
      integer JCAC(NNZAC)       ! work space
      double precision AAC(NNZAC) ! work space
      integer NNZFZZ            ! nnz in FZ per col. eqn.
      integer NNZFW             ! total nnz in FW
      integer NNZFWZ(NCOL)      ! nnz in FW for dz per col. eqn.
      integer NNZFWD(NCOL)      ! nnz in FW for dz per col. eqn.
                                !   (coming from del F/ del dz derivatives)
      integer NNZFWY(NCOL)      ! nnz in FW for y per col. eqn.
      integer NNZFWU(NCOL)      ! nnz in FW for u per col. eqn.
      integer NNZAW             ! nnz in AW
      integer NNZCW             ! total nnz in CW
      integer NNZCWW(NCOL)      ! nnz in CW for each col point
      integer NNZFU             ! total nnz in FU
      integer NNZFUZ(NCOL)      ! nnz in FU for dz per col. eqn.
      integer NNZFUD(NCOL)      ! nnz in FU for dz per col. eqn.
                                !   (coming from del F/ del dz derivatives)
      integer NNZFUY(NCOL)      ! nnz in FU for y per col. eqn.
      integer NNZFUU(NCOL)      ! nnz in FU for u per col. eqn.
      integer NNZFPP            ! nnz in FP per col. eqn.
      integer NNZAU             ! nnz in AU
      integer NNZCU             ! total nnz in CU
      integer NNZCUU(NCOL)      ! nnz in CU for each col point
      double precision WDUMMY(NCOL*(NZ+NY)) ! dummy
      double precision UDUMMY(NCOL*NU)      ! dummy

      integer i, k, l, j
      integer lcn, lcnoff, nioff, njoff, lw, ivaroff

C
C     Shouldn't be necessary, but initialize to zero
C
      do i = 1, 2*NZ+NY+NU+NP
         call DCOPY(NZ+NY, 0d0, 0, DF(1,i), 1)
      enddo
C
C     Get structure from model
C
      call DAEMODEL_DF_STRUC(NZ, NY, NU, NP, DF, LDDF)
C
C     Get structure for additional elemental constraints
C
      if( NAC.gt.0 ) then
         call ADDCON_DF(NZ, NY, NU, NCOL, NAC, WDUMMY, UDUMMY, NNZAC,
     1        IRAC, JCAC, AAC)
      endif
C
C     First equations for C_i
C


C
C     FZ (in order of collocation equations)
C
      lcn = 0
      do i = 1, NZ+NY
         do j = 1, NZ
            if( DF(i,j).ne.0d0 ) then
               lcn     = lcn + 1
               IR(lcn) = i
               JC(lcn) = j
            endif
         enddo
      enddo
      NNZFZZ = lcn
      nioff  = NZ+NY
      lcnoff = lcn - NNZFZZ
      do i = 2, NCOL
         do j = 1, NNZFZZ
            IR(lcn+j) = IR(lcnoff+j) + nioff
            JC(lcn+j) = JC(lcnoff+j)
         enddo
         lcn    = lcn    + NNZFZZ
         lcnoff = lcnoff + NNZFZZ
      enddo
C
C     FW (outer order: collocation equations)
C
      NNZFW = lcn
      do l = 1, NCOL
         nioff = (l-1)*(NZ+NY)
C
C     First all the stuff for the dz's
C
         NNZFWZ(l) = lcn
         njoff     = 0
         do k = 1, NCOL
            ivaroff = (k-1)*(NZ+NY)
            do j = 1, NZ
               lw = IVARC1(ivaroff+j)
               if( lw.le.NCOL*(NZ+NY)+NAC ) then
                  do i = 1, NY+NZ
                     if( DF(i,j).ne.0d0 ) then
                        lcn     = lcn + 1
                        IR(lcn) = i   + nioff
                        JC(lcn) = lw  + njoff
                     endif
                  enddo
               endif
            enddo
         enddo
         NNZFWZ(l) = lcn - NNZFWZ(l)
C
C     Now the stuff for the dz's for those  del F/del dz  derivatives
C
         NNZFWD(l) = lcn
         njoff     = 0
         ivaroff   = (l-1)*(NZ+NY)
         do j = 1, NZ
            lw = IVARC1(ivaroff+j)
            if( lw.le.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,NZ+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFWD(l) = lcn - NNZFWD(l)
C
C     Now the y's
C
         NNZFWY(l) = lcn
         njoff     = 0
         ivaroff   = NZ + (l-1)*(NZ+NY)
         do j = 1, NY
            lw = IVARC1(ivaroff+j)
            if( lw.le.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,2*NZ+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFWY(l) = lcn - NNZFWY(l)
C
C     Now the u's
C
         NNZFWU(l) = lcn
         ivaroff = NCOL*(NZ+NY) + (l-1)*NU
         do j = 1, NU
            lw = IVARC1(ivaroff+j)
            if( lw.le.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,2*NZ+NY+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFWU(l) = lcn - NNZFWU(l)
      enddo
      NNZFW = lcn - NNZFW
C
C     AW
C
      NNZAW = lcn
      nioff = NCOL*(NZ+NY)
      njoff = 0
      do k = 1, NNZAC
         j = JCAC(k)
         lw = IVARC1(j)
         if( lw.le.NCOL*(NZ+NY)+NAC ) then
            i       = IRAC(k)
            lcn     = lcn + 1
            IR(lcn) = i   + nioff
            JC(lcn) = lw  + njoff
         endif
      enddo
      NNZAW = lcn - NNZAW
C
C     CW
C
      NNZCW = lcn
      nioff = 0
      njoff = 0
      do k = 1, NCOL
         NNZCWW(k) = lcn
         ivaroff = (k-1)*(NZ+NY)
         do i = 1, NZ
            lw = IVARC1(ivaroff+i)
            if( lw.le.NCOL*(NZ+NY)+NAC ) then
               lcn     = lcn + 1
               IR(lcn) = i   + nioff
               JC(lcn) = lw  + njoff
            endif
         enddo
         NNZCWW(k) = lcn - NNZCWW(k)
      enddo
      NNZCW = lcn - NNZCW
C
C     And now N
C

C
C     FU (outer order: collocation equations)
C
      NNZFU = lcn
      do l = 1, NCOL
         nioff = (l-1)*(NZ+NY)
C
C     First all the stuff for the dz's (identities later...)
C
         NNZFUZ(l) = lcn
         njoff     = -NCOL*(NZ+NY)-NAC
         do k = 1, NCOL
            ivaroff = (k-1)*(NZ+NY)
            do j = 1, NZ
               lw = IVARC1(ivaroff+j)
               if( lw.gt.NCOL*(NZ+NY)+NAC ) then
                  do i = 1, NY+NZ
                     if( DF(i,j).ne.0d0 ) then
                        lcn     = lcn + 1
                        IR(lcn) = i   + nioff
                        JC(lcn) = lw  + njoff
                     endif
                  enddo
               endif
            enddo
         enddo
         NNZFUZ(l) = lcn - NNZFUZ(l)
C
C     Now the stuff for the dz's for those  del F/del dz  derivatives
C
         NNZFUD(l) = lcn
         njoff     = -NCOL*(NZ+NY)-NAC
         ivaroff   = (l-1)*(NZ+NY)
         do j = 1, NZ
            lw = IVARC1(ivaroff+j)
            if( lw.gt.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,NZ+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFUD(l) = lcn - NNZFUD(l)
C
C     Now the y's
C
         NNZFUY(l) = lcn
         njoff     = -NCOL*(NZ+NY)-NAC
         ivaroff   = (l-1)*(NZ+NY)+NZ
         do j = 1, NY
            lw = IVARC1(ivaroff+j)
            if( lw.gt.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,2*NZ+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFUY(l) = lcn - NNZFUY(l)
C
C     Now the u's
C
         NNZFUU(l) = lcn
         ivaroff   = NCOL*(NZ+NY) + (l-1)*NU
         do j = 1, NU
            lw = IVARC1(ivaroff+j)
            if( lw.gt.NCOL*(NZ+NY)+NAC ) then
               do i = 1, NY+NZ
                  if( DF(i,2*NZ+NY+j).ne.0d0 ) then
                     lcn     = lcn + 1
                     IR(lcn) = i   + nioff
                     JC(lcn) = lw  + njoff
                  endif
               enddo
            endif
         enddo
         NNZFUU(l) = lcn - NNZFUU(l)
      enddo
      NNZFU = lcn - NNZFU
C
C     AU
C
      NNZAU = lcn
      nioff = NCOL*(NZ+NY)
      njoff = -NCOL*(NZ+NY)-NAC
      do k = 1, NNZAC
         j = JCAC(k)
         lw = IVARC1(j)
         if( lw.gt.NCOL*(NZ+NY)+NAC ) then
            i       = IRAC(k)
            lcn     = lcn + 1
            IR(lcn) = i   + nioff
            JC(lcn) = lw  + njoff
         endif
      enddo
      NNZAU = lcn - NNZAU
C
C     FP
C
      NNZFPP = lcn
      njoff  = 0
      do i = 1, NZ+NY
         do j = 1, NP
            if( DF(i,2*NZ+NY+NU+j).ne.0d0 ) then
               lcn     = lcn + 1
               IR(lcn) = i
               JC(lcn) = j + njoff
            endif
         enddo
      enddo
      NNZFPP = lcn - NNZFPP
      nioff  = NZ+NY
      lcnoff = lcn - NNZFPP
      do i = 2, NCOL
         do j = 1, NNZFPP
            IR(lcn+j) = IR(lcnoff+j) + nioff
            JC(lcn+j) = JC(lcnoff+j)
         enddo
         lcn    = lcn    + NNZFPP
         lcnoff = lcnoff + NNZFPP
      enddo
C
C     CU
C
      NNZCU = lcn
      nioff = 0
      njoff = -NCOL*(NZ+NY)-NAC
      do k = 1, NCOL
         NNZCUU(k) = lcn
         ivaroff = (k-1)*(NZ+NY)
         do i = 1, NZ
            lw = IVARC1(ivaroff+i)
            if( lw.gt.NCOL*(NZ+NY)+NAC ) then
               lcn     = lcn + 1
               IR(lcn) = i   + nioff
               JC(lcn) = lw  + njoff
            endif
         enddo
         NNZCUU(k) = lcn - NNZCUU(k)
      enddo
      NNZCU = lcn - NNZCU

      if( lcn.ne.NNZA+NCOL*NZ ) then
c         write(*,*) 'lcn = ', lcn,' NNZA+NCOL*NZ = ',NNZA+NCOL*NZ
         stop
      endif
C
C     That's it
C
 9999 continue
      return
      end
