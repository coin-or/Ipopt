C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine PARTELE(NZ, NY, NU, NP, NCOL, NAC, NELE, NNZA, NNZAC,
     1     TI, X, IELE, IVARC, LRW, RW, LIW, IW, IERR)
C
C     $Id: partele.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Find elemental partitioning based on element IELE
C
C     Author:  Andreas Waechter
C              c/o Group of Larry Biegler
C              Department of Chemical Engineering
C              Carnegie Mellon University
C              Pittsburgh, PA
C
      implicit none

      integer NZ, NY, NU, NP, NCOL, NELE, NAC, IELE, NNZA, NNZAC
      double precision TI(NELE+1)
C      double precision X(NZ + NELE*(NCOL*(NZ+NY+NU)+NZ) + NP)
      double precision X(*)
      integer IVARC(NCOL*(NZ+NY+NU))
      integer LRW               ! at least max(20*NNZA + NCOL*(NZ+NY+NU),
                                !              NNZA + (NZ+NY)*(2*NZ+NY+NU+NP) +
                                !              NZ + NNZAC + NCOL*(NZ+NY+NU) )
      double precision RW(LRW+1)
      integer LIW               ! at least max(40*NNZA + 13*NCOL*(NZ+NY+NU),
                                !              40*NNZA + 10*NCOL + 2*NNZAC)
      integer IW(LIW+1)
      integer IERR

      integer p_nnzfwz, p_nnzfwy, p_nnzfwu, p_nnzcww, p_nnzfuz, p_nnzfuy
      integer p_nnzfuu, p_nnzcuu, nnzfzz, nnzfw, nnzcw, nnzfu, nnzfpp
      integer p_rwend, p_iwend, m, n, nnz, la, lddf, p_a, p_df, p_zval
      integer p_irn, p_jcn, lw, lu, lp, lc, ln, p_w, p_iw, k, i, l
      integer p_ikeep, nind, j, iflag, nnzcu, lan, p_cscale
      integer p_aac, p_wdummy, p_udummy, p_nnzfwd, p_nnzfud, p_irac
      integer p_jcac, nnzaw, nnzau
      double precision u
C
C     MA28 common block...
C
      double precision EPS, RMIN, RESID
      integer IRNCP, ICNCP, MINIRN, MINICN, IRANK, LP28, MP28
      logical ABORT1, ABORT2, LBLOCK, GROW
      common /MA28ED/ LP28, MP28, LBLOCK, GROW
      common /MA28FD/ EPS, RMIN, RESID, IRNCP, ICNCP, MINIRN, MINICN,
     1     IRANK, ABORT1, ABORT2

      IERR    = 0
      p_rwend = 0
      p_iwend = 0

      m    = NCOL*(NZ+NY) + NAC
      n    = NCOL*(NZ+NY+NU)
      la   = 20*NNZA            ! That should be PLENTY!
      lddf = NZ+NY
C
C     get some work space
C
      p_a      = p_rwend
      p_df     = p_a      + NNZA
      p_zval   = p_df     + lddf*(2*NZ+NY+NU+NP)
      p_aac    = p_zval   + NZ
      p_wdummy = p_aac    + NNZAC
      p_udummy = p_wdummy + NCOL*(NZ+NY)
      p_rwend  = p_udummy + NCOL*NU
      if( p_rwend.gt.LRW ) goto 9098

      p_irn    = p_iwend
      p_jcn    = p_irn    + la
      p_nnzfwz = p_jcn    + la
      p_nnzfwd = p_nnzfwz + NCOL
      p_nnzfwy = p_nnzfwd + NCOL
      p_nnzfwu = p_nnzfwy + NCOL
      p_nnzcww = p_nnzfwu + NCOL
      p_nnzfuz = p_nnzcww + NCOL
      p_nnzfud = p_nnzfuz + NCOL
      p_nnzfuy = p_nnzfud + NCOL
      p_nnzfuu = p_nnzfuy + NCOL
      p_nnzcuu = p_nnzfuu + NCOL
      p_irac   = p_nnzcuu + NCOL
      p_jcac   = p_irac   + NNZAC
      p_iwend  = p_jcac   + NNZAC
      if( p_iwend.gt.LIW ) goto 9099
C
C     First determine structure of jacobian (unpartitioned)
C
      do i = 1, n
         IVARC(i) = i
      enddo
      call STRUCTURE(NZ, NY, NU, NP, NCOL, NAC, IVARC, NNZA,
     1     IW(p_irn+1), IW(p_jcn+1), lddf, RW(p_df+1), NNZAC,
     2     IW(p_irac+1), IW(p_jcac+1), RW(p_aac+1),
     3     nnzfzz, nnzfw, IW(p_nnzfwz+1), IW(p_nnzfwd+1),
     4     IW(p_nnzfwy+1), IW(p_nnzfwu+1), nnzaw, nnzcw, IW(p_nnzcww+1),
     5     nnzfu, IW(p_nnzfuz+1), IW(p_nnzfud+1), IW(p_nnzfuy+1),
     6     IW(p_nnzfuu+1), nnzfpp, nnzau, nnzcu, IW(p_nnzcuu+1),
     7     RW(p_wdummy+1), RW(p_udummy+1))

      nnz  = nnzfw + nnzfu + NNZAC
C
C     First compute Jacobian of element IELE
C
      lw  = 1 + (IELE-1)*(NCOL*(NZ+NY) + NZ)
      lu  = 1 + NZ + NELE*(NCOL*(NZ+NY)+NZ) + (IELE-1)*NU
      lp  = 1 + NZ + NELE*(NCOL*(NZ+NY+NU)+NZ)
      lc  = 1
      ln  = 1 + NCOL*NNZFZZ + NNZFW + NNZAW + NNZCW
      lan = 1 + NCOL*NNZFZZ + NNZFW + NNZAW
      call FILLBLOK(NZ, NY, NU, NP, NCOL, NAC, NELE, X(lw), X(lu),
     1     X(lp), TI, IELE, IVARC, IVARC, nnzfzz, nnzfw, IW(p_nnzfwz+1),
     2     IW(p_nnzfwd+1), IW(p_nnzfwy+1), IW(p_nnzfwu+1), nnzfu,
     3     IW(p_nnzfuz+1), IW(p_nnzfud+1), IW(p_nnzfuy+1),
     4     IW(p_nnzfuu+1), nnzfpp, nnzaw, nnzau, IW(p_irn+lc),
     5     IW(p_jcn+lc), RW(p_a+lc), IW(p_irn+ln), IW(p_jcn+ln),
     6     RW(p_a+lan), lddf, RW(p_df+1), RW(p_zval+1), NNZAC,
     7     IW(p_irac+1), IW(p_jcac+1), RW(p_aac+1))
      p_rwend = p_a + la
C
C     Filter out unnecessary stuff to obtain the matrix that we wish
C     to partition
C
      do i = 1, NNZFW + NNZAW
         RW(p_a  +i) = RW(p_a  +NCOL*NNZFZZ+i)
         IW(p_irn+i) = IW(p_irn+NCOL*NNZFZZ+i)
         IW(p_jcn+i) = IW(p_jcn+NCOL*NNZFZZ+i)
      enddo
      ln  = ln  - 1
      lan = lan - 1
      do i = 1, NNZFU + NNZAU
         RW(p_a  +NNZFW+NNZAW+i) = RW(p_a  +lan+i)
         IW(p_irn+NNZFW+NNZAW+i) = IW(p_irn+ln +i)
         IW(p_jcn+NNZFW+NNZAW+i) = IW(p_jcn+ln +i) + m
      enddo
      p_rwend = p_a     + la
      p_iwend = p_jcn   + la
C
C     Scale the columns so that y's stay in basis if possible
C
      p_cscale = p_rwend
      p_rwend  = p_cscale + n
      if( p_rwend.gt.LRW ) goto 9098
      call DCOPY(n, 1d0, 0, RW(p_cscale+1), 1)
      k = NZ + 1
      do i = 1, NCOL
         call DCOPY(NY, 1d4, 0, RW(p_cscale+k), 1)
         k = k + NZ + NY
      enddo
      do i = 1, nnz
         RW(p_a+i) = RW(p_a+i)*RW(p_cscale+IW(p_jcn+i))
      enddo
      p_rwend = p_cscale
C
C     call MA28 to compute partition
C
      p_w     = p_rwend
      p_rwend = p_w     + n
      if( p_rwend.gt.LRW ) goto 9098

      p_ikeep = p_iwend
      p_iw    = p_ikeep + 5*n
      p_iwend = p_iw    + 8*n
      if( p_iwend.gt.LIW ) goto 9099

      ABORT1 = .false.
      ABORT2 = .false.
      LBLOCK = .false.
      MP28 = 0
      u = 1d0
      call MA28AD(n, nnz, RW(p_a+1), la, IW(p_irn+1), la, IW(p_jcn+1),
     1     u, IW(p_ikeep+1), IW(p_iw+1), RW(p_w+1), iflag)
      if( iflag.lt.0 .and. iflag.ne.-14 ) goto 9028
C
C     Get the partitioning out of IKEEP (old Claudia Schmid trick)
C
      k = 0
      do i = 1, n
         l = IW(p_ikeep+2*n+i)
         if( l.lt.0 ) then
C           indepentent variable
            k = k + 1
            IW(p_ikeep+k) = -l
         endif
      enddo
      nind = n-m
      if( k.ne.nind ) goto 9006
      k = m
      l = 0
      do i = 1, N
         do j = 1, nind
            if( i.eq.IW(p_ikeep+j) ) then
               k = k + 1
               IVARC(k) = i
               goto 110
            endif
         enddo
         l = l + 1
         IVARC(l) = i
 110     continue
      enddo
C
C     that's it
C
      goto 9999
C
C     Error handling
C
 9006 IERR = 6
      goto 9999
 9028 IERR = 28
      goto 9999
 9098 IERR = 98
      goto 9999
 9099 IERR = 99
      goto 9999
C
 9999 continue
      return
      end
