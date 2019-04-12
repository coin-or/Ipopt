C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine FILLBLOK(NZ, NY, NU, NP, NCOL, NAC, NELE, W, U, P, TI,
     1     IELE, IVARC, IVARC1, NNZFZZ, NNZFW, NNZFWZ, NNZFWD, NNZFWY,
     2     NNZFWU, NNZFU, NNZFUZ, NNZFUD, NNZFUY, NNZFUU, NNZFPP, NNZAW,
     3     NNZAU, IRC, JCC, C, IRN, JCN, N,
     3     LDDF, DF, ZVAL, NNZAC, IRAC, JCAC, AAC)
C
C     $Id: fillblok.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Compute nonzeros for one finite element
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none
      integer NZ, NY, NU, NP, NCOL, NAC, NELE
      double precision W(NZ + NCOL*(NZ+NY))
      double precision U(*)                  !(NCOL*NU)
      double precision P(NP)
      double precision TI(NELE+1)
      integer IELE
      integer IVARC(NCOL*(NZ+NY+NU))
      integer IVARC1(NCOL*(NZ+NY+NU)) ! inverse of IVARC (only needed if NAC>0)
      integer NNZFZZ            ! nnz in FZ per col. eqn.
      integer NNZFW             ! total nnz in FW
      integer NNZFWZ(NCOL)      ! nnz in FW for dz per col. eqn.
      integer NNZFWD(NCOL)      ! nnz in FW for dz (del F/del dz) per col. eqn.
      integer NNZFWY(NCOL)      ! nnz in FW for y per col. eqn.
      integer NNZFWU(NCOL)      ! nnz in FW for u per col. eqn.
      integer NNZFU             ! total nnz in FU
      integer NNZFUZ(NCOL)      ! nnz in FU for dz per col. eqn.
      integer NNZFUD(NCOL)      ! nnz in FU for dz (del F/del dz) per col. eqn.
      integer NNZFUY(NCOL)      ! nnz in FU for y per col. eqn.
      integer NNZFUU(NCOL)      ! nnz in FU for u per col. eqn.
      integer NNZFPP            ! nnz in FP per col. eqn.
      integer NNZAW
      integer NNZAU
      integer IRC(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer JCC(NCOL*NNZFZZ+NNZFW+NNZAW)
      double precision C(NCOL*NNZFZZ+NNZFW+NNZAW)
      integer IRN(NNZFU+NCOL*NNZFPP+NNZAU)
      integer JCN(NNZFU+NCOL*NNZFPP+NNZAU)
      double precision N(NNZFU+NCOL*NNZFPP+NNZAU)
      integer LDDF
      double precision DF(LDDF, 2*NZ+NY+NU+NP) ! work array
      double precision ZVAL(NZ) ! work array
      integer NNZAC
      integer IRAC(NNZAC)
      integer JCAC(NNZAC)
      double precision AAC(NNZAC)

      integer i, k, l, lfz, lfw, lfu, lfp
      integer nioff, lw, ivaroff, lwoff, jc, ir
      double precision h, hrho, trho, homega, dummy

      include 'DAE2NLP.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DAENLP/

      lfz = 0
      lfw = lfz + NNZFZZ*NCOL
      lfu = 0
      lfp = lfu + NNZFU + NNZAU
      h   = TI(IELE+1) - TI(IELE)
C
C     Loop over the collocation equations
C
      do l = 1, NCOL
         nioff = -(l-1)*(NZ+NY)
C
C     First compute values of z at current collocation point
C
         hrho = h * RHO(l)
         trho = TI(IELE) + hrho
         call APPROX(1, NZ, NY, NU, NELE, NCOL, COEF, W, U, TI, IELE,
     1        trho, OMEGAQ(1,l), dummy, ZVAL, dummy, dummy)
C
C     Now get derivatives from model
C
C     I don't know if we should make sure every time that DF is indeed
C     initialized to zero...  actually, i don't think it is necessary
C         call DCOPY(NZ, 0d0, 0, DF, LDDF+1)
         call DAEMODEL_DF(NZ, NY, NU, NP, trho, ZVAL,
     1        W(NZ+(l-1)*(NZ+NY)+1), W(NZ+(l-1)*(NZ+NY)+NZ+1),
     2        U((l-1)*NU+1), P, DF, LDDF)
C
C     FZ
C
         do i = 1, NNZFZZ
            C(lfz+i) = DF(IRC(lfz+i)+nioff,JCC(lfz+i))
         enddo
         lfz = lfz + NNZFZZ
C
C     FW
C
C     start with dz's
C
         ivaroff = 0
         lwoff   = 0
         k       = 1
         homega  = hrho*OMEGAQ(k,l)
         do i = 1, NNZFWZ(l)
            jc = JCC(lfw+i)
            ir = IRC(lfw+i) + nioff
 100        lw = IVARC(jc+ivaroff) + lwoff
            if( lw.gt.NZ ) then
               lwoff  = lwoff - NZ-NY
               k      = k + 1
               homega = hrho*OMEGAQ(k,l)
               goto 100
            endif
            C(lfw+i) = homega*DF(ir,lw)
         enddo
         lfw = lfw + NNZFWZ(l)
C
C     now  del F/del dz  entries for dz's
C
         ivaroff = 0
         lwoff   = -(l-1)*(NZ+NY)+NZ
         do i = 1, NNZFWD(l)
            jc = JCC(lfw+i)
            ir = IRC(lfw+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            C(lfw+i) = DF(ir,lw)
         enddo
         lfw = lfw + NNZFWD(l)
C
C     now y's
C
         ivaroff = 0
         lwoff   = -(l-1)*(NZ+NY)+NZ
         do i = 1, NNZFWY(l)
            jc = JCC(lfw+i)
            ir = IRC(lfw+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            C(lfw+i) = DF(ir,lw)
         enddo
         lfw = lfw + NNZFWY(l)
C
C     and the u's
C
         ivaroff = 0
         lwoff   = -NCOL*(NZ+NY) - (l-1)*NU + 2*NZ + NY
         do i = 1, NNZFWU(l)
            jc = JCC(lfw+i)
            ir = IRC(lfw+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            C(lfw+i) = DF(ir,lw)
         enddo
         lfw = lfw + NNZFWU(l)
C
C     FU
C
C     start with dz's
C
         ivaroff = NCOL*(NZ+NY)+NAC
         lwoff   = 0
         k       = 1
         homega  = hrho*OMEGAQ(k,l)
         do i = 1, NNZFUZ(l)
            jc = JCN(lfu+i)
            ir = IRN(lfu+i) + nioff
 110        lw = IVARC(jc+ivaroff) + lwoff
            if( lw.gt.NZ ) then
               lwoff  = lwoff - NZ-NY
               k      = k + 1
               homega = hrho*OMEGAQ(k,l)
               goto 110
            endif
            N(lfu+i) = homega*DF(ir,lw)
         enddo
         lfu = lfu + NNZFUZ(l)
C
C     now  del F/del dz  entries for dz's
C
         ivaroff = NCOL*(NZ+NY)+NAC
         lwoff   = -(l-1)*(NZ+NY) + NZ
         do i = 1, NNZFUD(l)
            jc = JCN(lfu+i)
            ir = IRN(lfu+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            N(lfu+i) = DF(ir,lw)
         enddo
         lfu = lfu + NNZFUD(l)
C
C     now y's
C
         ivaroff = NCOL*(NZ+NY)+NAC
         lwoff   = -(l-1)*(NZ+NY) + NZ
         do i = 1, NNZFUY(l)
            jc = JCN(lfu+i)
            ir = IRN(lfu+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            N(lfu+i) = DF(ir,lw)
         enddo
         lfu = lfu + NNZFUY(l)
C
C     and the u's
C
         ivaroff = NCOL*(NZ+NY)+NAC
         lwoff   = -NCOL*(NZ+NY) - (l-1)*NU + 2*NZ + NY
         do i = 1, NNZFUU(l)
            jc = JCN(lfu+i)
            ir = IRN(lfu+i) + nioff
            lw = IVARC(jc+ivaroff) + lwoff
            N(lfu+i) = DF(ir,lw)
         enddo
         lfu = lfu + NNZFUU(l)
C
C     FP
C
         do i = 1, NNZFPP
            jc = JCN(lfp+i)
            ir = IRN(lfp+i) + nioff
            N(lfp+i) = DF(ir,2*NZ+NY+NU+jc)
         enddo
         lfp = lfp + NNZFPP
C
      enddo
C
C     Nonzeros for additional constraints
C
      if( NAC.gt.0 ) then
         call ADDCON_DF(NZ, NY, NU, NCOL, NAC, W(NZ+1), U, NNZAC, IRAC,
     1        JCAC, AAC)
         do k = 1, NNZAC
            jc = JCAC(k)
            lw = IVARC1(jc)
            if( lw.le.NCOL*(NZ+NY)+NAC ) then
               lfw = lfw + 1
               C(lfw) = AAC(k)
            else
               lfu = lfu + 1
               N(lfu) = AAC(k)
            endif
         enddo
      endif

      return
      end
