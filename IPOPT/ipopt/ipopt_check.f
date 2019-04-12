C*******************************************************************************
C
      subroutine IPOPT_CHECK(NORIG, XORIG, M, NLBO, ILBO,
     1     BNDS_LO, NUBO, IUBO, BNDS_UO, V_LO, V_UO, LAM,
     1     LRW, RW, LIW, IW, ITER, IERR, EV_F, EV_C, EV_G, EV_A,
     1     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: ipopt_check.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Substitute for IPOPT to check derivatives by finite differences
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB    This subroutine is meant to allow an easy check of the derivatives
CB    provided by the user.  Since it has the same arguments are the main
CB    subroutine 'IPOPT', the check can be done simply by renaming the called
CB    routine.
CB    The derivatives are check around the given starting point (in XORIG)
CB    by using a simple first order finite difference test.  The size of the
CB    perturmation is given in the parameter 'DPIVTOL', and a possible error
CB    is reported, if the finite difference estimate deviates from the user
CB    provided derivatives value by more than 'DTOL' (in a relative measure).
CB    First derivatives are check in any case.  The Hessians are only check
CB    if IFULL = 1 and IQUASI = 0.
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter        10/03/02
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
CP   NORIG     I    INT    total number of variables
CP   XORIG     I    DP     point around with finite difference should be
CP                           performed
CP   M         I    INT    number of equality constraints
CP   NLBO      -    INT    (ignored)
CP   ILBO      -    INT    (ignored)
CP   BNDS_LO   -    DP     (ignored)
CP   NUBO      -    INT    (ignored)
CP   IUBO      -    INT    (ignored)
CP   BNDS_UO   -    DP     (ignored)
CP   V_LO      -    DP     (ignored)
CP   V_UO      -    DP     (ignored)
CP   LAM       W    DP     used as workspace
CP   LRW       I    INT    length of RW
CP   RW        W    DP     DP work space
CP   LIW       I    INT    length of IW
CP   IW        W    INT    INT work space
CP   IERR     I/O   INT    as in IPOPT
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
CS    DCOPY
CS    DAXPY
CS    CONSTR
CS    C_OUT
CS    GET_F
CS    GET_G
CS    GET_C
CS    GET_H
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
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NORIG
      double precision XORIG(NORIG)
      integer M
      integer NLBO
      integer ILBO(NLBO)
      double precision BNDS_LO(NLBO)
      integer NUBO
      integer IUBO(NUBO)
      double precision BNDS_UO(NUBO)
      double precision V_LO(NLBO)
      double precision V_UO(NUBO)
      double precision LAM(M)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer ITER, IERR	
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
      integer lrs_constr, lis_constr, lrw_constr, liw_constr, idummy
      integer kconstr(6), lrs_end, lis_end
      integer p_rwend, p_iwend, p_x, p_x2, p_c, p_c2, p_g, p_a
      integer p_cscale, p_ivar, p_avar, p_acon, p_hess, p_hess0
      integer p_g2, p_a2, p_irnh, p_icnh

      double precision delta, tol, h, f2, f, approx, reldiff, aa, dummy
      integer nz, i, j, k, l, nnzh
      integer errcount
      character*1 c_str
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
CTODO make work with dynamic memory allocation
CTODO adapt CARGS argument list
      errcount = 0

      delta = QPIVTOL
      tol   = QTOL
      iter  = 0

      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,' CHECKING DERIVATIVES:')
      call C_OUT(2,0,1,' ')
      write(line,*) 'Perturbation parameter (QPIVTOL) = ',delta
      call C_OUT(2,0,1,line)
      write(line,*) 'Error tolerance        (QTOL)    = ',tol
      call C_OUT(2,0,1,line)
C
C     Initialize CONSTR
C
      lrs_end = 0
      lis_end = 0
      if( M.ne.0 ) then
         lrw_constr = LRW
         liw_constr = LIW
         call CONSTR(0, idummy, NORIG, NORIG-M, M, idummy, idummy,
     1        idummy, NORIG, XORIG, dummy, dummy, dummy, nz,
     2        idummy, lrs_constr, dummy, lis_constr, idummy,
     3        lrw_constr, RW, liw_constr, IW, IERR, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'ipopt_check: CONSTR returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         kconstr(1) = lrs_constr
         kconstr(3) = lis_constr
         kconstr(5) = lrw_constr
         kconstr(6) = liw_constr

         kconstr(2) = lrs_end
         lrs_end    = kconstr(2) + lrs_constr
         kconstr(4) = lis_end
         lis_end    = kconstr(4) + lis_constr
      else
         nz = 0
      endif
      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,'----- Begin Check of First Derivatives  -----')
      write(line,1000) nz
 1000 format('Number of nonzeros in Jacobian: ',i8)
      call C_OUT(2,0,1,line)

      p_rwend = lrs_end
      p_iwend = lis_end

      p_x     = p_rwend
      p_x2    = p_x + NORIG
      p_c     = p_x2 + NORIG
      p_c2    = p_c  + M
      p_g     = p_c2 + M
      p_a     = p_g  + NORIG
      p_cscale= p_a  + nz
      if( QSCALE.ge.3 ) then
         p_rwend = p_cscale + M + NORIG
      elseif( QSCALE.ge.2 ) then
         p_rwend = p_cscale + M
      else
         p_rwend = p_cscale + 1
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      p_ivar  = p_iwend
      p_acon  = p_ivar + NORIG
      p_avar  = p_acon + nz
      p_iwend = p_avar + nz
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      call DCOPY(NORIG, XORIG, 1, RW(p_x+1), 1)
C
C     Assume no scaling
C
      if( QSCALE.ge.3 ) then
         call DCOPY(M+NORIG, 1.d0, 0, RW(p_cscale+1), 1)
      elseif( QSCALE.ge.2 ) then
         call DCOPY(M, 1.d0, 0, RW(p_cscale+1), 1)
      else
         RW(p_cscale+1) = 1.d0
      endif
C
C     No permutation of variables for derivative test
C
      do i = 1, NORIG
         IW(p_ivar+i) = i
      enddo
C
C     Compute values of objective function and constraints as well as their
C     first derivatives at given point
C
      call GET_F(NORIG, RW(p_x+1), IW(p_ivar+1), NORIG, XORIG,
     1     M, RW(p_cscale+1), f, EV_F, DAT, IDAT)
CUPDATE GET_F
         IERR = 3475
         goto 9999
      call GET_G(NORIG, RW(p_x+1), IW(p_ivar+1), NORIG, XORIG,
     1     M, RW(p_cscale+1), RW(p_g+1), LRW-p_rwend, RW(p_rwend+1),
     2     IERR, EV_G, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'ipopt_check: GET_G(a) returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( M.gt.0 ) then
         call GET_C(iter, NORIG, NORIG-M, RW(p_x+1), IW(p_ivar+1),
     1        NORIG, XORIG, M, RW(p_cscale+1), RW(p_c+1),
     2        KCONSTR, lrs_end, RW, lis_end, IW, LRW-p_rwend,
     3        RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     4        DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'ipopt_check: Error in GET_C(a), IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call CONSTR(10, iter, NORIG, NORIG-M, M, IW(p_ivar+1), 0,
     1        idummy, NORIG, XORIG, RW(p_cscale+1), RW(p_rwend+1),
     2        RW(p_a+1), IW(p_acon+1), IW(p_avar+1),
     3        kconstr(1), RW(kconstr(2)+1), kconstr(3),
     4        IW(kconstr(4)+1), LRW-p_rwend-2, RW(p_rwend+3),
     5        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'ipopt_check: Error in CONSTR(a), IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif

      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,' #var          : number of variable')
      call C_OUT(2,0,1,' #con          : number of constraint '//
     1'(0 for objective function)')
      call C_OUT(2,0,1,
     1' exact Deriv   : value given by derivative subroutine')
      call C_OUT(2,0,1,
     1' numeric Deriv : value obtained by finite difference')
      call C_OUT(2,0,1,
     1' rel. Diff     : relative difference')
      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,
     1'#var   #con    exact Deriv     numeric Deriv       rel. Diff')

      do i = 1, NORIG
         h = dmax1(delta, dabs(RW(p_x+i)*delta))
         call DCOPY(NORIG, RW(p_x+1), 1, RW(p_x2+1), 1)
         RW(p_x2+i) = RW(p_x2+i) + h
         call GET_F(NORIG, RW(p_x2+1), IW(p_ivar+1), NORIG, XORIG,
     1        M, RW(p_cscale+1), f2, EV_F, DAT, IDAT)
         call GET_C(iter, NORIG, NORIG-M, RW(p_x2+1), IW(p_ivar+1),
     1        NORIG, XORIG, M, RW(p_cscale+1), RW(p_c2+1),
     2        KCONSTR, lrs_end, RW, lis_end, IW, LRW-p_rwend,
     3        RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     4        DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'ipopt_check_derivs: Error in GET_C(b), IERR = ',IERR
            call C_OUT(2,0,1,line)
            IERR = 0
         endif
         approx = (f2-f)/h
         if( RW(p_g+i).ne.0.d0 ) then
            reldiff = abs((approx-RW(p_g+i))/RW(p_g+i))
         else
            reldiff = abs(approx)
         endif
         if( reldiff.gt.tol ) then
            c_str = '*'
            write(line,7000) i, 0, RW(p_g+i), approx, reldiff, c_str
            call C_OUT(2,0,1,line)
            errcount = errcount + 1
         else
            c_str = ' '
         endif
         do j = 1, M
            aa = 0.d0
            do k = 1, nz
               if( IW(p_acon+k).eq.j .and. IW(p_avar+k).eq.i ) then
                  aa = aa + RW(p_a+k)
               endif
            enddo
            approx = (RW(p_c2+j)-RW(p_c+j))/h
            if( aa.ne.0.d0 ) then
               reldiff = abs((approx-aa)/aa)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.tol ) then
               c_str = '*'
               write(line,7000) i, j, aa, approx, reldiff, c_str
               call C_OUT(2,0,1,line)
               errcount = errcount + 1
C            else
C               c_str = ' '
            endif
         enddo
      enddo

      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,'-----  End  Check of First Derivatives  -----')
      call C_OUT(2,0,1,' ')

      if( QFULL.ne.1 .or. QQUASI.ne.0 ) goto 6000
C
C     Check second derivatives
C
      call C_OUT(2,0,1,'----- Begin Check of Second Derivatives -----')
      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,' #con          : number of constraint '//
     1'(0 for objective function)')
      call C_OUT(2,0,1,' #var1         : number of first variable')
      call C_OUT(2,0,1,' #var2         : number of second variable')
      call C_OUT(2,0,1,
     1' exact Deriv   : value given by derivative subroutine')
      call C_OUT(2,0,1,
     1' numeric Deriv : value obtained by finite difference')
      call C_OUT(2,0,1,
     1' rel. Diff     : relative difference')
      call C_OUT(2,0,1,' ')
      call C_OUT(2,0,1,'#con   #var1 #var2   exact Deriv     '//
     1     'numeric Deriv       rel. Diff')
C     First get number of nonzeros in Hessian
      call GET_H(-1, NORIG, NORIG-M, 0, RW(p_x+1), IW(p_ivar+1),
     1     NORIG, XORIG, RW(p_cscale+1), M, nnzh, LAM, dummy,
     1     idummy, idummy,
     1     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1     IERR, EV_H, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'ipopt_check: GET_H(a) returns IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      call C_OUT(2,0,1,' ')
      write(line,1010) nnzh
 1010 format('Number of nonzeros in Hessian: ',i8)
      call C_OUT(2,0,1,line)
      call C_OUT(2,0,1,' ')
      endif

      p_hess  = p_rwend
      p_hess0 = p_hess  + nnzh
      p_g2    = p_hess0 + nnzh
      p_a2    = p_g2    + NORIG
      p_rwend = p_a2    + nz
      if( p_rwend.gt.LRW ) then
        IERR = 98
        goto 9999
      endif
      p_irnh  = p_iwend
      p_icnh  = p_irnh + nnzh
      p_iwend = p_icnh + nnzh
      if( p_iwend.gt.LIW ) then
        IERR = 99
        goto 9999
      endif

      do l = 0, M

         if( l.eq.0 ) then
            call C_OUT(2,1,1,
     1     '---- Begin Check of Hessian of Objective function -----')
         else
            write(line,1020) l
 1020       format('---- Begin Check of Hessien of Constraint',
     1           i6,' -----')
            call C_OUT(2,1,1,line)
         endif

         call DCOPY(M, 0.d0, 0, LAM, 1)
         if( l.gt.0 ) then
            LAM(l) = 1.d0
         endif
         call GET_H(iter, NORIG, NORIG-M, 0, RW(p_x+1), IW(p_ivar+1),
     1        NORIG, XORIG, RW(p_cscale+1), M, nnzh, LAM, RW(p_hess+1),
     1        IW(p_irnh+1), IW(p_icnh+1), LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_H, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'ipopt_check: GET_H(b) returns IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( l.eq.0 ) then
            call DCOPY(nnzh, RW(p_hess+1), 1, RW(p_hess0+1), 1)
         else
            call DAXPY(nnzh, -1.d0, RW(p_hess0+1), 1, RW(p_hess+1), 1)
         endif

         do i = 1, NORIG
            h = dmax1(delta, dabs(RW(p_x+i)*delta))
            call DCOPY(NORIG, RW(p_x+1), 1, XORIG, 1)
            XORIG(i) = XORIG(i) + h

            if( l.eq.0 ) then
               call DCOPY(NORIG, XORIG, 1, RW(p_x2+1), 1)
               call GET_G(NORIG, RW(p_x2+1), IW(p_ivar+1), NORIG, XORIG,
     1              M, RW(p_cscale+1), RW(p_g2+1), LRW-p_rwend,
     2              RW(p_rwend+1), IERR, EV_G, DAT, IDAT)
               if( IERR.ne.0 ) then
                  write(line,*)
     1                 'ipopt_check: GET_G(a) returns IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
            else
               call CONSTR(10, iter, NORIG, NORIG-M, M, IW(p_ivar+1), 0,
     1              idummy, NORIG, XORIG, RW(p_cscale+1), RW(p_rwend+1),
     2              RW(p_a2+1), IW(p_acon+1), IW(p_avar+1),
     3              kconstr(1), RW(kconstr(2)+1), kconstr(3),
     4              IW(kconstr(4)+1), LRW-p_rwend-2, RW(p_rwend+3),
     5              LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5              EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
               if( IERR.ne.0 ) then
                  write(line,*)
     1                 'ipopt_check: Error in CONSTR(b), IERR = ',IERR
                  call C_OUT(2,0,1,line)
                  goto 9999
               endif
               call DCOPY(NORIG, 0.d0, 0, RW(p_g +1), 1)
               call DCOPY(NORIG, 0.d0, 0, RW(p_g2+1), 1)
               do k = 1, nz
                  if( IW(p_acon+k).eq.l ) then
                     RW(p_g +IW(p_avar+k)) = RW(p_a +k)
                     RW(p_g2+IW(p_avar+k)) = RW(p_a2+k)
                  endif
               enddo
            endif

            do j = 1, NORIG
               aa = 0.d0
               do k = 1, nnzh
                  if( (IW(p_icnh+k).eq.i .and. IW(p_irnh+k).eq.j) .or.
     1                (IW(p_icnh+k).eq.j .and. IW(p_irnh+k).eq.i)
     1                 ) then
                     aa = aa + RW(p_hess+k)
                  endif
               enddo
               approx = (RW(p_g2+j)-RW(p_g+j))/h
               if( aa.ne.0.d0 ) then
                  reldiff = abs((approx-aa)/aa)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.tol ) then
                  c_str = '*'
                  write(line,7010) l, i, j, aa, approx, reldiff, c_str
                  call C_OUT(2,0,1,line)
                  errcount = errcount + 1
               endif
            enddo

         enddo

      enddo

 6000 continue
      call C_OUT(2,0,1,' ')
      if( errcount.eq.0 ) then
         write(line,1030) tol
 1030    format('All derivatives are within specified',
     1        ' specified tolerance of', d9.2,'.')
         call C_OUT(2,0,1,line)
      else
         write(line,1040) errcount, tol
 1040    format('In',i6,' cases the derivatives are not within',
     1        ' specified tolerance of', d9.2,'.')
         call C_OUT(2,0,1,line)
      endif

 9999 continue
      return

 7000 format(i5,i6,d17.8,d17.8,d17.8,a2)
 7010 format(i5,i6,i6,d17.8,d17.8,d17.8,a2)
      end
