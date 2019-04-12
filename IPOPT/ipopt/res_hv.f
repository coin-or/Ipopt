C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine RES_HV(ITER, N, NIND, M, NORIG, XORIG, CSCALE, IVAR,
     1     NFIX, IFIX, SOFLAG, C, VIN, VOUT,
     1     NLB, ILB, NUB, IUB, S_L, S_U,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: res_hv.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute product VOUT = Hessian  times  VIN
CT    Here, Hessian = A A^T    if SOFLAG = .false. ,and
CT    Hessian = A A^T + sum( c_i  \na_xx c_i ) otherwise
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      implicit none
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer ITER
      integer N
      integer NIND
      integer M
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      logical SOFLAG
      double precision C(M)
      double precision VIN(N)
      double precision VOUT(N)
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision S_L(NLB)
      double precision S_U(NUB)
      integer KCONSTR(6)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, p_tmp, p_vino, p_vouto, p_xtmp
      integer i, idummy, k
      double precision vnrm, h, DNRM2, v
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_rwend = 0
      p_iwend = 0
      IERR = 0
C
C     First compute  A A^T v
C

      p_tmp = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif

      call CONSTR(9, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, VIN, RW(p_tmp+1), idummy, idummy,
     3     KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'res_hv: Warning in CONSTR-9, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'res_hv: Error in CONSTR-9, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

      call CONSTR(8, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, RW(p_tmp+1), VOUT, idummy, idummy,
     3     KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4     IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'res_hv: Warning in CONSTR-8, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'res_hv: Error in CONSTR-8, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

      p_rwend = p_tmp
C
C     If wanted, also add the term involving the Hessian of c^T c
C
      if( SOFLAG ) then

         if( QSCALE.eq.2 .or. CSCALE(1).ne.1.d0 ) then
            write(line,*) 'res_hv: can''t work with scaling yet'
            call C_OUT(2,0,1,line)
            IERR = 4
            goto 9999
         endif

         p_vino  = p_rwend
         p_vouto = p_vino + NORIG
         p_rwend = p_vouto + NORIG
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         if( NFIX.ne.0 ) then
            call DCOPY(NORIG, 0.d0, 0, RW(p_vino+1), 1)
         endif
         do i = 1, N
            RW(p_vino+IVAR(i)) = VIN(i)
         enddo
C
C     The following IF block is copied from get_hv (unify?)
C
         if( QHESS.eq.2 ) then
C
C     Estimate Hessian vector product by finite difference
C     NOTE:  We assume that this will be premultiplied by Z, so that we don't
C     comute A * VIN  here!!!!!
C
C     Determine step size
C
c               i    = IDAMAX(NORIG, RW(p_vino+1), 1)
c               vnrm = dabs(RW(p_vino+i))
            vnrm = DNRM2(NORIG, RW(p_vino+1), 1)
            if( vnrm.eq.0.d0 ) then ! VOUT is already set to 0
               goto 9999
            else
               h = 1d-6/vnrm
               h = 1d-6/dmin1(vnrm, 1d0)
c               h = 1d-6
C               h = 1d-1
C
C     Make sure that no bounds are violated
C
c                  goto 1444
               do i = 1, NLB
                  k = ILB(i)
                  v = RW(p_vino+IVAR(k))
                  if( v.lt.0.d0 ) then
                     h = dmin1(h, -S_L(i)/v)
                  endif
               enddo
               do i = 1, NUB
                  k = IUB(i)
                  v = RW(p_vino+IVAR(k))
                  if( v.gt.0.d0 ) then
                     h = dmin1(h, S_U(i)/v)
                  endif
               enddo
 1444          continue
               write(*,*) 'vnrm = ',vnrm,', h = ',h
C
C     Compute vouto = A(x+h VIN) * C
C
               p_xtmp  = p_rwend
               p_rwend = p_xtmp + NORIG
               if( p_rwend.gt.LRW ) then
                  IERR  = 98
                  goto 9999
               endif
               call DCOPY(NORIG, XORIG, 1, RW(p_xtmp+1), 1)
               call DAXPY(NORIG, h, RW(p_vino+1), 1, RW(p_xtmp+1), 1)
C
               call CONSTR(13, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1              NORIG, RW(p_tmp+1), CSCALE, C, RW(p_vouto+1),
     3              idummy, idummy, KCONSTR(1), RS(KCONSTR(2)+1),
     4              KCONSTR(3), IS(KCONSTR(4)+1), LRW-p_rwend,
     5              RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2              EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     6              EV_HCV, DAT, IDAT)
               if( IERR.lt.0 ) then
                  write(line,*)
     1                 'res_hv: Warning in CONSTR-13, IERR = ',IERR
                  call C_OUT(2,0,1,line)
               elseif( IERR.ne.0 ) then
                  write(line,*)
     1                 'res_hv: Error in CONSTR-13, IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
               p_rwend = p_xtmp
               call DSCAL(NORIG, 1/h, RW(p_vouto+1), 1)
            endif

         else
            call CONSTR(20, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, C, RW(p_vino+1), RW(p_vouto+1),
     3           idummy, idummy, KCONSTR(1), RS(KCONSTR(2)+1),
     4           KCONSTR(3), IS(KCONSTR(4)+1), LRW-p_rwend,
     5           RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2           EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     6           EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*)
     1              'res_hv: Warning in CONSTR-20, IERR = ',IERR
               call C_OUT(2,0,1,line)
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'res_hv: Error in CONSTR-20, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
         endif
         do i = 1, N
            VOUT(i) = VOUT(i) + RW(p_vouto+IVAR(i))
         enddo

      endif

 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine RES_HV_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)
      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)

      call CONSTR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      if( QHESS.eq.2 ) then
         LRW = LRW + max(M,2*N)
      else
         LRW = LRW + max(M,3*N)
      endif

      return
      end

