C Copyright (C) 2003, International Business Machines and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_SCALE(ITER, NORIG, XORIG, N, X, NIND, M, IVAR,
     1     NFIX, IFIX, NLBO, ILBO, BNDS_LO, NUBO, IUBO, BNDS_UO, CSCALE,
     1     KCONSTR, LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR, EV_F,
     1     EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: get_scale.F 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute scaling factors
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CBC
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      10/21/03  First version
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   ITER      I    INT    ignored, unless -1:
CP                           NIND is NNZA (nonzeros in Jacobian)
CP                           NFIX is NNZH (nonzeros in Hessian)
CP   NORIG     I    INT    total number of variables (incl. fixed vars)
CP   XORIG     I    INT    original starting point (unscaled)
CP   N         I    INT    number of nonfixed variables
CP   X         W    DP     work space
CP   NIND      I    INT    number of degrees of freedom
CP   M         I    INT    number of equality constraints
CP   IVAR      W    INT    work space
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   NLBO      I    INT    number of lower bounds (incl. fixed vars)
CP   ILBO      I    INT    indices of lower bounds (incl. fixed vars)
CP                            (e.g. BNDS_LO(i) is bound for XORIG(ILBO(i)) )
CP   BNDS_LO   I    DP     values of lower bounds (incl. fixed vars)
CP   NUBO      I    INT    number of upper bounds (incl. fixed vars)
CP   IUBO      I    INT    indices of upper bounds (incl. fixed vars)
CP                            (e.g. BNDS_UO(i) is bound for XORIG(IUBO(i)) )
CP   BNDS_UO   I    DP     values of upper bounds (incl. fixed vars)
CP   CSCALE    O    DB     scaling factors:
CP                           first M entries are scaling factors d_c for
CP                             constraints (c_alg = d_c * c_orig)
CP                           last N entries are scaling factors d_x for
CP                             variables (x_orig = d_x * x_alg)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS        S    DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS        S    INT    INT storage space (all!)
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective fuction
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
      include 'IPOPT.INC'

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer ITER
      integer NORIG
      double precision XORIG(NORIG)
      integer N
      double precision X(N)
      integer NIND
      integer M
      integer IVAR(N)
      integer NFIX
      integer IFIX(NFIX)
      integer NLBO
      integer ILBO(NLBO)
      double precision BNDS_LO(NLBO)
      integer NUBO
      integer IUBO(NUBO)
      double precision BNDS_UO(NUBO)
      double precision CSCALE(*)
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
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, p_cg, p_g, p_a, p_acon, p_avar, p_r, p_c
      integer i, idummy, ix, lfix, nz, nza, p_w, ifail, p_tmp, p_scale
      integer p_ilbo1, p_iubo1, retval, nv
      logical fixed, new
      double precision fscale_orig, basval, dummy(2), cscalmin, scalval
      double precision pl, pu, bl, bu, diff
      character*1 cdummy

      integer IDAMAX, FFINITE, READ_ENTRY

      character*80 line(3)

      integer NNZA, NNZH
      save    NNZA, NNZH
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( ITER.eq.-1 ) then
         NNZA = NIND
         NNZH = NFIX
         goto 9999
      endif
         
      IERR = 0
      p_iwend = 0
      p_rwend = 0

      if( MEMDBG ) then
         write(line,1)'get_scale', LRW, LIW
 1       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      if( QSCALE.eq.0 ) then
         CSCALE(1) = 1.d0
         goto 9999
      endif

      if( QSCALE.eq.-1 ) then
         CSCALE(1) = QSCALECUT
         goto 9999
      endif

      if( QSCALE.eq.5 ) then
C     Read scaling factors from file
         open(80,file='SCALE.DAT',status='old',err=9008)
         new = .true.
         do i = 1, N+M
            retval = READ_ENTRY(80, new, 'd', CSCALE(i), idummy, cdummy)
            if( retval.ne.0 ) goto 9009
         enddo
         close(80)
         if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
            write(line,8003) QFSCALE
            call C_OUT(1,0,2,line)
            write(line,8001)
            call C_OUT(1,0,3,line)
            do i = 1, M
               write(line,8002) i,CSCALE(i)
               call C_OUT(1,0,1,line)
            enddo
            write(line,8004)
            call C_OUT(1,0,3,line)
            do i = 1, N
               write(line,8005) i, CSCALE(M+i)
               call C_OUT(1,0,1,line)
            enddo
         endif
         goto 9999
      endif
C
C     Create dummy IVAR and X (moved into bounds)
C
      ix = 0
      lfix = 1
      p_ilbo1 = p_iwend
      p_iubo1 = p_ilbo1 + NORIG
      p_iwend = p_iubo1 + NORIG
      do i = 1, NORIG
         IW(p_ilbo1+i) = 0
         IW(p_iubo1+i) = 0
      enddo
      do i = 1, NLBO
         IW(p_ilbo1+ILBO(i)) = i
      enddo
      do i = 1, NUBO
         IW(p_iubo1+IUBO(i)) = i
      enddo
      do i = 1, NORIG
         fixed = .false.
         if( lfix.le.NFIX ) then
            if( IFIX(lfix).eq.i ) then
               lfix = lfix + 1
               fixed = .true.
            endif
         endif
         if( .not.fixed ) then
            ix = ix + 1
            IVAR(ix) = i
            pl = 0.d0
            pu = 0.d0
            if( IW(p_ilbo1+i).ne.0 ) then
               bl = BNDS_LO(IW(p_ilbo1+i))
               pl = QTHETA*max(1.d0,abs(bl))
            endif
            if( IW(p_iubo1+i).ne.0 ) then
               bu = BNDS_UO(IW(p_iubo1+i))
               pu = QTHETA*max(1.d0,abs(bu))
            endif
            X(ix) = XORIG(i)
            if( pl.ne.0.d0 ) then
               if( pu.ne.0.d0 ) then
                  diff = BNDS_UO(IW(p_iubo1+i)) - BNDS_LO(IW(p_ilbo1+i))
                  X(ix) = dmax1(dmin1(bl + QFRAC*diff, bl + pl),
     1                 dmin1(dmax1(bu - QFRAC*diff, bu - pu), X(ix)))
               else
                  X(ix) = dmax1(bl + pl, X(ix))
               endif
            elseif( pu.ne.0.d0 ) then
               X(ix) = dmin1(bu - pu, X(ix))
            endif
         endif
      enddo
      p_iwend = p_ilbo1

      if( QSCALE.eq.1 .or. QSCALE.eq.2 ) then
C
C     Scale objective and constraints only
C
         fscale_orig = QFSCALE
         QFSCALE = 1d0
         p_g     = p_rwend
         p_rwend = p_g + N
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     Compute gradient of objective and use max element to compute scaling for F
C
         call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, 0, idummy,
     1        dummy, 0, idummy, dummy, 0.d0, RW(p_g+1),
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1        IERR, EV_G, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'get_scale: get_g returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         i = IDAMAX(N, RW(p_g+1), 1)
         QFSCALE = abs(RW(p_g+i))
         write(line,*) 'get_scale: |g|_inf = ',QFSCALE
         call C_OUT(2,2,1,line)
         if( FFINITE(QFSCALE).eq.0 ) then
            call C_OUT(2,0,1, 'Objective function gradient cannot be eva
     1luated at starting point!')
            IERR = 11
            goto 9999
         endif
         QFSCALE = QSCALECUT/dmax1(QSCALECUT,QFSCALE)
         QFSCALE = QFSCALE*fscale_orig
         write(line,*) 'get_scale: QFSCALE = ',QFSCALE
         call C_OUT(2,2,1,line)
         p_rwend = p_g
C
C     If there are constraints, compute their gradients
C
         if( M.gt.0 ) then
            if( QSCALE.eq.1 ) then
               CSCALE(1) = 1.d0
            else
               call DCOPY(M, 1.d0, 0, CSCALE, 1)
            endif
C     Need to choose ITER = -100 so that A will be recomputed next time
C     ( CONSTR stores A in a scaled manner)
            p_cg    = p_rwend
            p_rwend = p_cg + M
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            call CONSTR(14, -100, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, RW(p_cg+1),
     2           idummy, idummy,
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.lt.0 ) then
               write(line,*)
     1              'get_scale: Warning in CONSTR-14, IERR = ', IERR
               call C_OUT(2,0,1,line)
            elseif( IERR.ne.0 ) then
               write(line,*)
     1              'get_scale: Error in CONSTR-14, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            basval = QSCALECUT
            if( QSCALE.eq.1 ) then
               i = IDAMAX(M, RW(p_cg+1), 1)
               CSCALE(1) = abs(RW(p_cg+i))
               write(line,*) 'get_scale: |A|_inf = ',CSCALE(1)
               call C_OUT(2,2,1,line)
               CSCALE(1) = basval/dmax1(basval, CSCALE(1))
               write(line,*) 'get_scale: CSCALE  = ',CSCALE(1)
               call C_OUT(2,2,1,line)
            else
               cscalmin = 1.d300
               do i = 1, M
                  CSCALE(i) = basval/dmax1(basval, RW(p_cg+i))
                  cscalmin = min(cscalmin, CSCALE(i))
               enddo
               write(line,*) 'get_scale: smallest CSCALE = ',cscalmin
               call C_OUT(2,2,1,line)
               if( cscalmin.eq.1.d0 ) then
                  call C_OUT(2,2,1,
     1                 'get_scale: No scaling of constraints necessary')
                  QSCALE = 1
               endif
               if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
                  write(line,8001)
 8001             format(/,
     1                'Information about constraint scaling factors:',/)
                  call C_OUT(1,0,3,line)
                  do i = 1, M
                     write(line,8002) i,CSCALE(i)
 8002                format('Scaling for constraint ',i8,': ',d10.4)
                     call C_OUT(1,0,1,line)
                  enddo
               endif
            endif
            p_rwend = p_cg
         endif
      elseif( .false. .and. (QSCALE.eq.3 .or. QSCALE.eq.4) ) then
C
C     Compute equilibration factors
C
         call DCOPY(M+N, 1.d0, 0, CSCALE, 1)
         fscale_orig = QFSCALE
         basval = abs(fscale_orig)
         if( QSCALE.eq.3 ) then
            scalval = basval
         else
            scalval = 1.d0
         endif
         QFSCALE = 1.d0
         nza = NNZA
C
C     Compute (unscaled) Jacobian of the constraints
C     (Workspace should be large enough if memory allocation works, otherwise
C     we could have a problem here).
C
         p_a     = p_rwend
         p_rwend = p_a + nza + N
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_acon  = p_iwend
         p_avar  = p_acon + nza + N
         p_iwend = p_avar + nza + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         if( M.gt.0 ) then
            call CONSTR(10, -100, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, RW(p_a+1),
     2           IW(p_acon+1), IW(p_avar+1),
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if(  IERR.ne.0 ) then
               write(line,*)
     1              'get_scale: Error in CONSTR-10, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            NNZA = int(dummy(1))
C
C     Take out those values above the threshold value and rescale
C
            nz = 0
            do i = 1, nza
               if( abs(RW(p_a+i)).gt.basval ) then
                  nz = nz + 1
                  RW(p_a+nz)    = RW(p_a+i)/scalval
                  IW(p_acon+nz) = IW(p_acon+i)
                  IW(p_avar+nz) = IW(p_avar+i)
               endif
            enddo
         else
            nz = 0
         endif
C
C     Now add the gradient for the objective function
C
         p_rwend = p_a + nz + N
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         nza = nz
         call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, 0, idummy,
     1        dummy, 0, idummy, dummy, 0.d0, RW(p_a+nza+1),
     1        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1        IERR, EV_G, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'get_scale: get_g returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         do i = 1, N
            if( abs(RW(p_a+nza+i)).gt.basval ) then
               nz = nz + 1
               RW(p_a+nz)    = RW(p_a+nza+i)/scalval
               IW(p_acon+nz) = M+1
               IW(p_avar+nz) = i
            endif
         enddo
         p_rwend = p_a + nz
         if( nz.eq.0 ) then
            if( fscale_orig.lt.0.d0 ) then
               QFSCALE = -1.d0
            else
               QFSCALE = 1.d0
            endif
            call DCOPY(M, 1.d0, 0, CSCALE, 1)
            call C_OUT(2,2,1,'get_scale: No scaling necessary.')
         else
            p_r     = p_rwend
            p_c     = p_r + M+1
            p_w     = p_c + N
            p_rwend = p_w + 2*(M+1) + 3*N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
#ifdef HAVE_MC29
C
C     It seems that the second version scales the constraints usually not up
C
C            call MC29AD(M+1, N, nz, RW(p_a+1), IW(p_acon+1),
C     1           IW(p_avar+1), RW(p_r+1), RW(p_c+1), RW(p_w+1), QCNR,
C     2           ifail)
            call MC29AD(N, M+1, nz, RW(p_a+1), IW(p_avar+1),
     1           IW(p_acon+1), RW(p_c+1), RW(p_r+1), RW(p_w+1), QCNR,
     2           ifail)
            if( ifail.ne.0 ) then
               write(line,*) 'get_scale: MC29AD returned IFAIL = ',ifail
               call C_OUT(2,0,1,line)
               IERR = 416
               goto 9999
            endif
#else
            call C_OUT(2,0,1,
     1       'Error: This scaling option is only available with MC29.')
            call C_OUT(2,0,1,
     2'       However, this version of IPOPT is compiled without MC29')
            IERR = 4
            goto 9999
#endif
            do i = 1, M
               CSCALE(i) = exp(RW(p_r+i))
            enddo
            QFSCALE = exp(RW(p_r+M+1))
            if( fscale_orig.lt.0.d0 ) QFSCALE = -QFSCALE
            do i = 1, N
               CSCALE(M+i) = exp(RW(p_c+i))
            enddo
            if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
               write(line,8003) QFSCALE
 8003          format(/'Scaling factor for objective function = ',d10.4)
               call C_OUT(1,0,2,line)
               write(line,8001)
               call C_OUT(1,0,3,line)
               do i = 1, M
                  write(line,8002) i,CSCALE(i)
                  call C_OUT(1,0,1,line)
               enddo
               write(line,8004)
 8004          format(/,
     1              'Information about variable scaling factors:',/)
               call C_OUT(1,0,3,line)
               do i = 1, N
                  write(line,8005) i, CSCALE(M+i)
 8005             format('Scaling for variable ',i8,': ',d10.4)
                  call C_OUT(1,0,1,line)
               enddo
            endif
         endif   
      elseif( .false. .and.(QSCALE.eq.5 .or. QSCALE.eq.6) ) then
         if( QFULL.eq.0 .or. QQUASI.ne.0 ) then
            call C_OUT(2,0,1,'get_scale: need exact derivatives and full
     1 space for this scaling option.')
            IERR = 4
            goto 9999
         endif
C
C     Compute equilibration factors
C
         call DCOPY(M+N, 1.d0, 0, CSCALE, 1)
         fscale_orig = QFSCALE
         basval = abs(fscale_orig)
         if( QSCALE.eq.5 ) then
            scalval = basval
         else
            scalval = 1.d0
         endif
         QFSCALE = 1.d0
C
C     First get Hessian of objective function
C
         p_a     = p_rwend
         p_tmp   = p_a + NNZH
         p_rwend = p_tmp + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         p_acon  = p_iwend
         p_avar  = p_acon + NNZA + NNZH
         p_iwend = p_avar + NNZA + NNZH
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call DCOPY(M, 0.d0, 0, RW(p_tmp+1), 1)
         call GET_H(-100, N, NIND, NFIX, X, IVAR, NORIG, XORIG,
     1        CSCALE, M, NNZH, RW(p_tmp+1), RW(p_a+1), IW(p_acon+1),
     1        IW(p_avar+1), LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
CDELETEME
         call DCOPY(NNZH, 0.d0, 0, RW(p_a+1), 1)
         p_rwend = p_a + NNZH + NNZA
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     Compute (unscaled) Jacobian of the constraints
C
         if( M.gt.0 ) then
            call CONSTR(10, -100, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, RW(p_a+NNZH+1),
     2           IW(p_acon+NNZH+1), IW(p_avar+NNZH+1),
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if(  IERR.ne.0 ) then
               write(line,*)
     1              'get_scale: Error in CONSTR-10, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
         endif
         NNZA = int(dummy(1))
C
C     Take out those values above the threshold value and rescale
C
         nz = 0
         do i = 1, NNZH
            if( abs(RW(p_a+i)).gt.basval ) then
               nz = nz + 1
               RW(p_a+nz)    = RW(p_a+i)/scalval
               IW(p_acon+nz) = IW(p_acon+i)
               IW(p_avar+nz) = IW(p_avar+i)
            endif
         enddo
         do i = NNZH+1, NNZH+NNZA
            if( abs(RW(p_a+i)).gt.basval ) then
               nz = nz + 1
               RW(p_a+nz)    = RW(p_a+i)/scalval
               IW(p_acon+nz) = IW(p_acon+i) + N
               IW(p_avar+nz) = IW(p_avar+i)
            endif
         enddo
         p_rwend = p_a + nz
         p_iwend = p_avar + nz
         if( nz.eq.0 ) then
            if( fscale_orig.lt.0.d0 ) then
               QFSCALE = -1.d0
            else
               QFSCALE = 1.d0
            endif
            call DCOPY(M, 1.d0, 0, CSCALE, 1)
            call C_OUT(2,2,1,'get_scale: No scaling necessary.')
         else
            p_scale = p_rwend
            p_rwend = p_scale + M+N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C
C     call MC30, MC29 or MC19 to obtain scaling parameters
C
#ifdef HAVE_MC30
            p_w     = p_rwend
            p_rwend = p_w     + 4*(N+M)
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__,
     1              'p_rwend = ',p_rwend
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            endif
            call MC30AD(N+M, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), RW(p_w+1), QCNR, ifail)
            p_rwend = p_w
            if( ifail.ne.0 ) then
               write(line,*) 'MC30AD returns ifail = ',ifail
               call C_OUT(2,0,1,line)
               IERR = 517
               goto 9999
            endif
#else
#ifdef HAVE_MC29
            call MC29_CALL(N+M, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), LRW-p_rwend,
     4           RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR)
            if( IERR.ne.0 ) then
               write(line,*) 'MC29call returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
#else
#ifdef HAVE_MC19
C
            p_w      = p_rwend  ! This will be REAL work space in MC19_CALL
            p_rwend  = p_w + (7*(N+M))/2 + 1
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            endif
            call MC19_CALL(N+M, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), RW(p_w+1),
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
            p_rwend = p_w
            if( IERR.ne.0 ) then
               write(line,*) 'MC19call returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
#else
            IERR = 97
            call C_OUT(2,0,1,
     1           'get_scale:  Compiled without MC30, MC29 or MC19.')
            call C_OUT(2,0,1,
     1           '            Choose a different value for ISCALE.')
            goto 9999
#endif
#endif
#endif
            do i = 1, M
               CSCALE(i) = exp(RW(p_scale+N+i))
            enddo
            if( fscale_orig.lt.0.d0 ) QFSCALE = -1.d0
            do i = 1, N
               CSCALE(M+i) = exp(RW(p_scale+i))
            enddo
            p_rwend = p_a
            if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
               write(line,8003) QFSCALE
               call C_OUT(1,0,2,line)
               write(line,8001)
               call C_OUT(1,0,3,line)
               do i = 1, M
                  write(line,8002) i,CSCALE(i)
                  call C_OUT(1,0,1,line)
               enddo
               write(line,8004)
               call C_OUT(1,0,3,line)
               do i = 1, N
                  write(line,8005) i, CSCALE(M+i)
                  call C_OUT(1,0,1,line)
               enddo
            endif
         endif   
      elseif( QSCALE.ge.3 .and. QSCALE.le.4 ) then
         if( (QFULL.eq.0 .or. QQUASI.ne.0) .and. QSCALE.eq.4 ) then
            call C_OUT(2,0,1,'get_scale: need exact derivatives and full
     1 space for this scaling option.')
            IERR = 4
            goto 9999
         endif
C
C     Compute equilibration factors
C
         call DCOPY(M+N, 1.d0, 0, CSCALE, 1)
         fscale_orig = QFSCALE
         basval = QSCALECUT
         scalval = 1.d0
         QFSCALE = 1.d0

         nz = 0
         if( QSCALE.eq.4 ) then
C
C     First get Hessian of objective function
C
            p_a     = p_rwend
            p_tmp   = p_a + NNZH
            p_rwend = p_tmp + M
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_acon  = p_iwend
            p_avar  = p_acon + NNZA + NNZH
            p_iwend = p_avar + NNZA + NNZH
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            call DCOPY(M, 0.d0, 0, RW(p_tmp+1), 1)
            call GET_H(-100, N, NIND, NFIX, X, IVAR, NORIG, XORIG,
     1           CSCALE, M, NNZH, RW(p_tmp+1), RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), LRW-p_rwend, RW(p_rwend+1),
     1           LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
            p_rwend = p_a + NNZH + NNZA
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            nz = NNZH
            nv = N+M
         else
            p_a     = p_rwend
            p_rwend = p_a + NNZA + N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
            p_acon  = p_iwend
            p_avar  = p_acon + NNZA + N
            p_iwend = p_avar + NNZA + N
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            nv = N+M+1
         endif
C
C     Compute (unscaled) Jacobian of the constraints
C
         if( M.gt.0 ) then
            call CONSTR(10, -100, N, NIND, M, IVAR, NFIX, IFIX,
     1           NORIG, XORIG, CSCALE, dummy, RW(p_a+nz+1),
     2           IW(p_acon+nz+1), IW(p_avar+nz+1),
     3           KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4           IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if(  IERR.ne.0 ) then
               write(line,*)
     1              'get_scale: Error in CONSTR-10, IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            NNZA = int(dummy(1))
         else
            NNZA = 0
         endif
C
C     Correct the position entries for Jacobian
C
         do i = 1, NNZA
            IW(p_acon+nz+i) = IW(p_acon+nz+i) + N
         enddo
         nz = nz + NNZA

         if( QSCALE.eq.3 ) then
C
C     Get gradient for objective function
C
            call GET_G(N, X, IVAR, NORIG, XORIG, M, CSCALE, 0, idummy,
     1           dummy, 0, idummy, dummy, 0.d0, RW(p_a+1),
     1           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1           IERR, EV_G, DAT, IDAT)
            if( IERR.ne.0 ) then
               write(line,*) 'get_scale: get_g returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            do i = 1, N
               IW(p_acon+nz+i) = M+N+1
               IW(p_avar+nz+i) = i
            enddo
            nz = nz + N
         endif
C
C     Take out those values above the threshold value and rescale
C
         nza = nz
         nz = 0
         do i = 1, nza
            if( abs(RW(p_a+i)).gt.basval ) then
               nz = nz + 1
               RW(p_a+nz)    = RW(p_a+i)/scalval
               IW(p_acon+nz) = IW(p_acon+i)
               IW(p_avar+nz) = IW(p_avar+i)
            endif
         enddo
         p_rwend = p_a + nz
         p_iwend = p_avar + nz
         if( nz.eq.0 ) then
            QFSCALE = fscale_orig
            call DCOPY(M, 1.d0, 0, CSCALE, 1)
            call C_OUT(2,2,1,'get_scale: No scaling necessary.')
         else
            p_scale = p_rwend
            p_rwend = p_scale + M+N
            if( p_rwend.gt.LRW ) then
               IERR = 98
               goto 9999
            endif
C
C     call MC30, MC29 or MC19 to obtain scaling parameters
C
#ifdef HAVE_MC30
            p_w     = p_rwend
            p_rwend = p_w     + 4*nv
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__,
     1              'p_rwend = ',p_rwend
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            endif
            call MC30AD(nv, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), RW(p_w+1), QCNR, ifail)
            p_rwend = p_w
            if( ifail.ne.0 ) then
               write(line,*) 'MC30AD returns ifail = ',ifail
               call C_OUT(2,0,1,line)
               IERR = 517
               goto 9999
            endif
#else
#ifdef HAVE_MC29
            call MC29_CALL(nv, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), LRW-p_rwend,
     4           RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR)
            if( IERR.ne.0 ) then
               write(line,*) 'MC29call returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
#else
#ifdef HAVE_MC19
C
            p_w      = p_rwend  ! This will be REAL work space in MC19_CALL
            p_rwend  = p_w + (7*(nv))/2 + 1
            if( p_rwend.gt.LRW ) then
#ifdef PRINT__LINE__
               write(line,*) 'Out of DP memory in line ',__LINE__
               call C_OUT(2,0,1,line)
#endif
               IERR = 98
               goto 9999
            endif
            call MC19_CALL(nv, nz, RW(p_a+1), IW(p_acon+1),
     1           IW(p_avar+1), RW(p_scale+1), RW(p_w+1),
     2           LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR)
            p_rwend = p_w
            if( IERR.ne.0 ) then
               write(line,*) 'MC19call returns IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
#else
            IERR = 97
            call C_OUT(2,0,1,
     1           'get_scale:  Compiled without MC30, MC29 or MC19.')
            call C_OUT(2,0,1,
     1           '            Choose a different value for ISCALE.')
            goto 9999
#endif
#endif
#endif
            do i = 1, M
               CSCALE(i) = exp(RW(p_scale+N+i))
            enddo
            if( fscale_orig.lt.0.d0 ) QFSCALE = -1.d0
            do i = 1, N
               CSCALE(M+i) = exp(RW(p_scale+i))
            enddo
            QFSCALE = fscale_orig
            if( QSCALE.eq.3 ) then
               QFSCALE = QFSCALE*exp(RW(p_scale+N+M+1))
               write(line,*) 'get_scale: QFSCALE = ',QFSCALE
               call C_OUT(2,2,1,line)
            endif
            p_rwend = p_a
            p_iwend = p_acon
            if( QCNR.gt.0 .and. QPRINT.ge.3 ) then
               write(line,8003) QFSCALE
               call C_OUT(1,0,2,line)
               write(line,8001)
               call C_OUT(1,0,3,line)
               do i = 1, M
                  write(line,8002) i,CSCALE(i)
                  call C_OUT(1,0,1,line)
               enddo
               write(line,8004)
               call C_OUT(1,0,3,line)
               do i = 1, N
                  write(line,8005) i, CSCALE(M+i)
                  call C_OUT(1,0,1,line)
               enddo
            endif
         endif   
      else
         call C_OUT(2,0,1,'initpoint: invalid value of QSCALE.')
         IERR = 4
         goto 9999
      endif
 9999 continue
      return

 9008 IERR = 8
      call C_OUT(2,0,1,'Error while trying to open SCALE.DAT')
      goto 9999
 9009 IERR = 8
      call C_OUT(2,0,1,'Error while reading SCALE.DAT')
      goto 9999
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine GET_SCALE_WS(N, M, NLB, NUB, NZA, NZH, LRW, LIW, DAT,
     1     IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, NZH, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer lrw1, liw1, nz
      character*80 line(1)

      LRW = 0
      LIW = 0

      if( QSCALE.eq.1 .or. QSCALE.eq.2 ) then
         call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
         LRW = max(LRW, lrw1+N)
         LIW = max(LIW, liw1)
         if( M.gt.0 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1+N)
            LIW = max(LIW, liw1)
         endif
      elseif( .false. .and. (QSCALE.eq.3 .or. QSCALE.eq.4) ) then
         call GET_G_WS(N, M, NLB, NUB, NZA, lrw1, liw1)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
         if( M.gt.0 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
         endif
         LRW = max(LRW,3*(M+1)+4*N) ! MC29
         LRW = LRW + (NZA+N)
         LIW = LIW + 2*(NZA+N)
      elseif( .false. .and. (QSCALE.eq.5 .or. QSCALE.eq.6) ) then
         call GET_H_WS(N, M, NLB, NUB, NZA, NZH, lrw1, liw1)
         LRW = max(LRW, lrw1)
         LIW = max(LIW, liw1)
         if( M.gt.0 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
         endif
#ifdef HAVE_MC30
         lrw1 = 4*(N+M)
#else
#ifdef HAVE_MC29
         call MC29_CALL_WS(N, M, NLB, NUB, NZA, NZH+NZA, lrw1, liw1)
#else
#ifdef HAVE_MC19
         call MC19_CALL_WS(N, M, NLB, NUB, NZA, NZH+NZA, lrw1, liw1)
         lrw1 = lrw1 + (7*(N+M))/2 + 1
#endif
#endif
#endif
         LRW = max(LRW, (N+M)+lrw1)
         LIW = max(LIW, liw1)

         LRW = LRW + (NZA+NZH)
         LIW = LIW + 2*(NZA+NZH)
      elseif( QSCALE.eq.3 .or. QSCALE.eq.4 ) then
         if( QSCALE.eq.4 ) then
            call GET_H_WS(N, M, NLB, NUB, NZA, NZH, lrw1, liw1)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
            nz = NZH
         else
            nz = N
         endif
         if( M.gt.0 ) then
            call CONSTR_WS(N, M, NLB, NUB, NZA, lrw1, liw1, DAT, IDAT)
            LRW = max(LRW, lrw1)
            LIW = max(LIW, liw1)
            nz = nz + NZA
         endif
#ifdef HAVE_MC30
         lrw1 = 4*(N+M+1)
#else
#ifdef HAVE_MC29
         call MC29_CALL_WS(N, M, NLB, NUB, NZA, nz, lrw1, liw1)
#else
#ifdef HAVE_MC19
         call MC19_CALL_WS(N, M, NLB, NUB, NZA, nz, lrw1, liw1)
         lrw1 = lrw1 + (7*(N+M+1))/2 + 1
#endif
#endif
#endif
         LRW = max(LRW, (N+M+1)+lrw1)
         LIW = max(LIW, liw1)

         LRW = LRW + (nz)
         LIW = LIW + 2*(nz)
      endif

      if( QSCALE.ge.1 ) then
         LIW = max(LIW, 2*N)
      endif

      if( QPRINT.ge.4 ) then
         write(line,1000)'get_scale_ws', LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end
