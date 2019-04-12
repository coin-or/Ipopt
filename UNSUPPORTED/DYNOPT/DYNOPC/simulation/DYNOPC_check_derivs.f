C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine FMAIN(MODNAM, STUBNAM)
C
C     Check Derivative Information from Model routines / DYNOPC version
C
C     $Id: DYNOPC_check_derivs.f 531 2004-03-11 01:31:07Z andreasw $
C
C     Author:    Andreas Waechter
C                Department of Chemical Engineering
C                Carnegie Mellon University
C                Pittsburgh, PA
C                USA
C
C     Date       09-30-01
C
      implicit none
      character*(*) MODNAM
      character*(*) STUBNAM
C
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/

      integer LDDF, LDRS
      parameter( LDDF = NZMAX+NYMAX )
      parameter( LDRS = NZMAX+NYMAX+NUMAX+NPMAX )

      double precision Z(NZMAX), ZB(2,NZMAX), DZ(NZMAX)
      double precision Y(NYMAX), YB(2,NYMAX)
      double precision U(NUMAX), UB(2,NUMAX)
      double precision P(NPMAX), PB(2,NPMAX)
      double precision F(NZMAX+NYMAX)

      double precision HZ(NZMAX), HDZ(NZMAX)
      double precision HY(NYMAX)
      double precision HU(NUMAX)
      double precision HP(NPMAX)
      double precision Z0(NZMAX), DZ0(NZMAX)
      double precision Y0(NYMAX)
      double precision U0(NUMAX)
      double precision P0(NPMAX)
      double precision F0(NZMAX+NYMAX)

      double precision DF (LDDF, 2*NZMAX+NYMAX+NUMAX+NPMAX)
      double precision DF0(LDDF, 2*NZMAX+NYMAX+NUMAX+NPMAX)
      double precision DF_STRUC(LDDF, 2*NZMAX+NYMAX+NUMAX+NPMAX)
      double precision LAM(NZMAX+NYMAX)
      double precision RHS(LDRS, 2*NZMAX+NYMAX+NUMAX+NPMAX)
      double precision RES(LDRS, 2*NZMAX+NYMAX+NUMAX+NPMAX)

      double precision VLARGE
      parameter( VLARGE = 1d20 )

      double precision EPS, EPS2, T, TOL
      parameter( EPS = 1.d-7, TOL = 1.d-6 )

      integer i, j, k, index, ierr, nlines
      character*80 filename
      character*3 c_str
      character*2 nz_str, name
      double precision exact, reldiff, approx
      double precision tol_daspk, tol_discret, tol_opt
      character*80 lines(7)

      EPS2 = dsqrt(EPS)
      T    = 0.d0
C
C     Read the *.cmp file
C
      call READ_CMP(MODNAM, STUBNAM, tol_daspk, tol_discret,
     1     tol_opt, zb, z, yb, y, ub, u, pb, p, ierr, nlines, lines)
      if( ierr.ne.0 ) then
         OutputUnit_COUT = 6
         call C_OUT(0,0,nlines,lines)
         stop
      endif

C
C     Obtain slightly perturbed starting point (to get rid of inital zeros)
C
      call PERTURB(NZ, Z, Z0, ZB, HZ, 'Z', EPS)
      call PERTURB(NY, Y, Y0, YB, HY, 'Y', EPS)
      call PERTURB(NU, U, U0, UB, HU, 'U', EPS)
      call PERTURB(NP, P, P0, PB, HP, 'P', EPS)
      do i = 1, NZ
         DZ0(i) = dble(i)*0.2d0
         HDZ(i) = EPS
      enddo
C
C     For ADOL-C version:  Write Tape
C
C     Initialize ADOL-C Model (i.e. create tape)
      call GET_FILENAME(MODNAM, MODNAM, ' ', filename)
!DEC$ ATTRIBUTES C, REFERENCE, ALIAS:'_model_init_' :: MODEL_INIT
      call MODEL_INIT(NZ, NY, NU, NP, TI(1), z0, dz, y0, u0, p0, f,
     1     filename)
C
C     Get Jacobian Structure
C
      do i = 1, 2*NZ+NY+NU+NP
         call DCOPY( NZ+NY, 0.d0, 0, DF_STRUC(1,i), 1 )
      enddo
      call MODEL_DF_STRUC(NZ, NY, NU, NP, DF_STRUC, LDDF)
C
C     Get Function and Jacobian at check point
C
      call MODEL_F(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0, F0)
      do i = 1, 2*NZ+NY+NU+NP
         call DCOPY( NZ+NY, 0.d0, 0, DF0(1,i), 1 )
      enddo
      call MODEL_DF(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0, DF0,
     1     LDDF)
C
C     Check Jacobian
C
      write(*,*)
      write(*,*)
     1' ==============================================================='
      write(*,*)
     1'                        CHECKING JACOBIAN '
      write(*,*)
     1' ==============================================================='

      call DCOPY(NZ, Z0, 1, Z, 1)
      do j = 1, NZ
C         write(*,*)
         Z(j) = Z0(j) + HZ(j)
         call MODEL_F(NZ, NY, NU, NP, T, Z, DZ0, Y0, U0, P0, F)
         do i = 1, NZ+NY
            exact  = DF0(i,j)
            approx = (F(i)-F0(i))/HZ(j)
            if( exact.ne.0.d0 ) then
               reldiff = abs((approx-exact)/exact)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.TOL ) then
               c_str = '<<<'
            else
               c_str = ' '
            endif
            if( DF_STRUC(i,j).ne.0.d0 ) then
               nz_str = 'nz'
            else
               nz_str = ' '
               if( exact.ne.0.d0 ) then
                  c_str = '<nz'
               endif
            endif
            if( c_str.ne.' ' )
     1      write(*,7000) 'Z',j,i, nz_str, exact, approx, reldiff, c_str
         enddo
         Z(j) = Z0(j)
      enddo
C      write(*,*)
C      write(*,*)
C     1' ---------------------------------------------------------------'
C      write(*,*)

      call DCOPY(NZ, DZ0, 1, DZ, 1)
      do j = 1, NZ
C         write(*,*)
         DZ(j) = DZ0(j) + HDZ(j)
         call MODEL_F(NZ, NY, NU, NP, T, Z0, DZ, Y0, U0, P0, F)
         do i = 1, NZ+NY
            exact  = DF0(i,NZ+j)
            approx = (F(i)-F0(i))/HDZ(j)
            if( exact.ne.0.d0 ) then
               reldiff = abs((approx-exact)/exact)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.TOL ) then
               c_str = '<<<'
            else
               c_str = ' '
            endif
            if( DF_STRUC(i,NZ+j).ne.0.d0 ) then
               nz_str = 'nz'
            else
               nz_str = ' '
               if( exact.ne.0.d0 ) then
                  c_str = '<nz'
               endif
            endif
            if( c_str.ne.' ' )
     1      write(*,7000) 'DZ',j,i, nz_str, exact, approx, reldiff,
     1           c_str
         enddo
         DZ(j) = DZ0(j)
      enddo
C      write(*,*)
C      write(*,*)
C     1' ---------------------------------------------------------------'
C      write(*,*)
      
      call DCOPY(NY, Y0, 1, Y, 1)
      do j = 1, NY
C         write(*,*)
         Y(j) = Y0(j) + HY(j)
         call MODEL_F(NZ, NY, NU, NP, T, Z0, DZ0, Y, U0, P0, F)
         do i = 1, NZ+NY
            exact  = DF0(i,2*NZ+j)
            approx = (F(i)-F0(i))/HY(j)
            if( exact.ne.0.d0 ) then
               reldiff = abs((approx-exact)/exact)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.TOL ) then
               c_str = '<<<'
            else
               c_str = ' '
            endif
            if( DF_STRUC(i,2*NZ+j).ne.0.d0 ) then
               nz_str = 'nz'
            else
               nz_str = ' '
               if( exact.ne.0.d0 ) then
                  c_str = '<nz'
               endif
            endif
            if( c_str.ne.' ' )
     1      write(*,7000) 'Y',j,i, nz_str, exact, approx, reldiff, c_str
         enddo
         Y(j) = Y0(j)
      enddo
C      write(*,*)
C      write(*,*)
C     1' ---------------------------------------------------------------'
C      write(*,*)

      call DCOPY(NU, U0, 1, U, 1)
      do j = 1, NU
C         write(*,*)
         U(j) = U0(j) + HU(j)
         call MODEL_F(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U, P0, F)
         do i = 1, NZ+NY
            exact  = DF0(i,2*NZ+NY+j)
            approx = (F(i)-F0(i))/HU(j)
            if( exact.ne.0.d0 ) then
               reldiff = abs((approx-exact)/exact)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.TOL ) then
               c_str = '<<<'
            else
               c_str = ' '
            endif
            if( DF_STRUC(i,2*NZ+NY+j).ne.0.d0 ) then
               nz_str = 'nz'
            else
               nz_str = ' '
               if( exact.ne.0.d0 ) then
                  c_str = '<nz'
               endif
            endif
            if( c_str.ne.' ' )
     1      write(*,7000) 'U',j,i, nz_str, exact, approx, reldiff, c_str
         enddo
         U(j) = U0(j)
      enddo
C      write(*,*)
C      write(*,*)
C     1' ---------------------------------------------------------------'
C      write(*,*)

      call DCOPY(NP, P0, 1, P, 1)
      do j = 1, NP
C         write(*,*)
         P(j) = P0(j) + HP(j)
         call MODEL_F(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P, F)
         do i = 1, NZ+NY
            exact  = DF0(i,2*NZ+NY+NU+j)
            approx = (F(i)-F0(i))/HP(j)
            if( exact.ne.0.d0 ) then
               reldiff = abs((approx-exact)/exact)
            else
               reldiff = abs(approx)
            endif
            if( reldiff.gt.TOL ) then
               c_str = '<<<'
            else
               c_str = ' '
            endif
            if( DF_STRUC(i,2*NZ+NY+NU+j).ne.0.d0 ) then
               nz_str = 'nz'
            else
               nz_str = ' '
               if( exact.ne.0.d0 ) then
                  c_str = '<nz'
               endif
            endif
            if( c_str.ne.' ' )
     1      write(*,7000) 'P',j,i, nz_str, exact, approx, reldiff, c_str
         enddo
         P(j) = P0(j)
      enddo
C
      write(*,*)
      write(*,*) 
     1     'Do you also want to check the Hessian vector products?'
      write(*,*)
      pause
C
C     Check Hessians
C
      call DCOPY(NZ+NY, 0.d0, 0, LAM, 1)
      do k = 1, NZ+NY
         LAM(k) = 3.d0
         write(*,*)
         write(*,*)
     1' ==============================================================='
         write(*,7001) k
 7001    format('                CHECKING HESSIAN FOR EQUATION ',i4)
         write(*,*)
     1' ==============================================================='
         write(*,*)
C
C     Compute Products of that Hessian for all Z entries
C
         if( NZ.le.0 ) goto 1500
         do j = 1, NZ
            call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RHS(1,j), 1)
            RHS(j,j) = 2.d0
         enddo
         call MODEL_HV( NZ, NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0,
     1        LAM, RHS, RES, LDRS)
         do j = 1, NZ
            call DSCAL(2*NZ+NY+NU+NP, 1.d0/6.d0, RES(1,j), 1)
         enddo

         do j = 1, NZ
            Z(j) = Z0(j) + HZ(j)
            call MODEL_DF(NZ, NY, NU, NP, T, Z, DZ0, Y0, U0, P0, DF,
     1           LDDF)
            do i = 1, 2*NZ+NY+NU+NP
               exact  = RES(i,j)
               approx = (DF(k,i)-DF0(k,i))/HZ(j)
               if( exact.ne.0.d0 ) then
                  reldiff = abs((approx-exact)/exact)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.TOL ) then
                  c_str = '<<<'
               else
                  c_str = ' '
               endif
               if( i.le.NZ ) then
                  name = 'Z'
                  index = i
               elseif( i.le.2*NZ ) then
                  name = 'DZ'
                  index = i - NZ
               elseif( i.le.2*NZ+NY ) then
                  name = 'Y'
                  index = i - 2*NZ
               elseif( i.le.2*NZ+NY+NU ) then
                  name = 'U'
                  index = i - 2*NZ-NY
               else
                  name = 'P'
                  index = i - 2*NZ-NY-NU
               endif
               if( c_str.ne.' ' )
     1         write(*,7010) k,name,index,'Z',j, exact, approx,
     1              reldiff, c_str
            enddo
C            write(*,*)
            Z(j) = Z0(j)
         enddo
C         write(*,*)
C     1' ---------------------------------------------------------------'
C         write(*,*)
C
C     Compute Products of that Hessian for all DZ entries
C
         do j = 1, NZ
            call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RHS(1,j), 1)
            RHS(NZ+j,j) = 2.d0
         enddo
         call MODEL_HV( NZ, NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0,
     1        LAM, RHS, RES, LDRS)
         do j = 1, NZ
            call DSCAL(2*NZ+NY+NU+NP, 1.d0/6.d0, RES(1,j), 1)
         enddo

         do j = 1, NZ
            DZ(j) = DZ0(j) + HDZ(j)
            call MODEL_DF(NZ, NY, NU, NP, T, Z0, DZ, Y0, U0, P0, DF,
     1           LDDF)
            do i = 1, 2*NZ+NY+NU+NP
               exact  = RES(i,j)
               approx = (DF(k,i)-DF0(k,i))/HDZ(j)
               if( exact.ne.0.d0 ) then
                  reldiff = abs((approx-exact)/exact)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.TOL ) then
                  c_str = '<<<'
               else
                  c_str = ' '
               endif
               if( i.le.NZ ) then
                  name = 'Z'
                  index = i
               elseif( i.le.2*NZ ) then
                  name = 'DZ'
                  index = i - NZ
               elseif( i.le.2*NZ+NY ) then
                  name = 'Y'
                  index = i - 2*NZ
               elseif( i.le.2*NZ+NY+NU ) then
                  name = 'U'
                  index = i - 2*NZ-NY
               else
                  name = 'P'
                  index = i - 2*NZ-NY-NU
               endif
               if( c_str.ne.' ' )
     1         write(*,7010) k,name,index,'DZ',j, exact, approx,
     1              reldiff, c_str
            enddo
C            write(*,*)
            DZ(j) = DZ0(j)
         enddo
C         write(*,*)
C     1' ---------------------------------------------------------------'
C         write(*,*)
C
C     Compute Products of that Hessian for all Y entries
C
 1500    if( NY.le.0 ) goto 1510
         do j = 1, NY
            call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RHS(1,j), 1)
            RHS(2*NZ+j,j) = 2.d0
         enddo
         call MODEL_HV( NY, NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0,
     1        LAM, RHS, RES, LDRS)
         do j = 1, NY
            call DSCAL(2*NZ+NY+NU+NP, 1.d0/6.d0, RES(1,j), 1)
         enddo

         do j = 1, NY
            Y(j) = Y0(j) + HY(j)
            call MODEL_DF(NZ, NY, NU, NP, T, Z0, DZ0, Y, U0, P0, DF,
     1           LDDF)
            do i = 1, 2*NZ+NY+NU+NP
               exact  = RES(i,j)
               approx = (DF(k,i)-DF0(k,i))/HY(j)
               if( exact.ne.0.d0 ) then
                  reldiff = abs((approx-exact)/exact)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.TOL ) then
                  c_str = '<<<'
               else
                  c_str = ' '
               endif
               if( i.le.NZ ) then
                  name = 'Z'
                  index = i
               elseif( i.le.2*NZ ) then
                  name = 'DZ'
                  index = i - NZ
               elseif( i.le.2*NZ+NY ) then
                  name = 'Y'
                  index = i - 2*NZ
               elseif( i.le.2*NZ+NY+NU ) then
                  name = 'U'
                  index = i - 2*NZ-NY
               else
                  name = 'P'
                  index = i - 2*NZ-NY-NU
               endif
               if( c_str.ne.' ' )
     1         write(*,7010) k,name,index,'Y',j, exact, approx,
     1              reldiff, c_str
            enddo
C            write(*,*)
            Y(j) = Y0(j)
         enddo
C         write(*,*)
C     1' ---------------------------------------------------------------'
C         write(*,*)
C
C     Compute Products of that Hessian for all U entries
C
 1510    if( NU.le.0 ) goto 1520
         do j = 1, NU
            call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RHS(1,j), 1)
            RHS(2*NZ+NY+j,j) = 2.d0
         enddo
         call MODEL_HV( NU, NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0,
     1        LAM, RHS, RES, LDRS)
         do j = 1, NU
            call DSCAL(2*NZ+NY+NU+NP, 1.d0/6.d0, RES(1,j), 1)
         enddo

         do j = 1, NU
            U(j) = U0(j) + HU(j)
            call MODEL_DF(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U, P0, DF,
     1           LDDF)
            do i = 1, 2*NZ+NY+NU+NP
               exact  = RES(i,j)
               approx = (DF(k,i)-DF0(k,i))/HU(j)
               if( exact.ne.0.d0 ) then
                  reldiff = abs((approx-exact)/exact)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.TOL ) then
                  c_str = '<<<'
               else
                  c_str = ' '
               endif
               if( i.le.NZ ) then
                  name = 'Z'
                  index = i
               elseif( i.le.2*NZ ) then
                  name = 'DZ'
                  index = i - NZ
               elseif( i.le.2*NZ+NY ) then
                  name = 'Y'
                  index = i - 2*NZ
               elseif( i.le.2*NZ+NY+NU ) then
                  name = 'U'
                  index = i - 2*NZ-NY
               else
                  name = 'P'
                  index = i - 2*NZ-NY-NU
               endif
               if( c_str.ne.' ' )
     1         write(*,7010) k,name,index,'U',j, exact, approx,
     1              reldiff, c_str
            enddo
C            write(*,*)
            U(j) = U0(j)
         enddo
C         write(*,*)
C     1' ---------------------------------------------------------------'
C         write(*,*)
C
C     Compute Products of that Hessian for all P entries
C
 1520    if( NP.le.0 ) goto 1530
         do j = 1, NP
            call DCOPY(2*NZ+NY+NU+NP, 0.d0, 0, RHS(1,j), 1)
            RHS(2*NZ+NY+NU+j,j) = 2.d0
         enddo
         call MODEL_HV( NP, NZ, NY, NU, NP, T, Z0, DZ0, Y0, U0, P0,
     1        LAM, RHS, RES, LDRS)
         do j = 1, NP
            call DSCAL(2*NZ+NY+NU+NP, 1.d0/6.d0, RES(1,j), 1)
         enddo

         do j = 1, NP
            P(j) = P0(j) + HP(j)
            call MODEL_DF(NZ, NY, NU, NP, T, Z0, DZ0, Y0, U, P0, DF,
     1           LDDF)
            do i = 1, 2*NZ+NY+NU+NP
               exact  = RES(i,j)
               approx = (DF(k,i)-DF0(k,i))/HP(j)
               if( exact.ne.0.d0 ) then
                  reldiff = abs((approx-exact)/exact)
               else
                  reldiff = abs(approx)
               endif
               if( reldiff.gt.TOL ) then
                  c_str = '<<<'
               else
                  c_str = ' '
               endif
               if( i.le.NZ ) then
                  name = 'Z'
                  index = i
               elseif( i.le.2*NZ ) then
                  name = 'DZ'
                  index = i - NZ
               elseif( i.le.2*NZ+NY ) then
                  name = 'Y'
                  index = i - 2*NZ
               elseif( i.le.2*NZ+NY+NU ) then
                  name = 'U'
                  index = i - 2*NZ-NY
               else
                  name = 'P'
                  index = i - 2*NZ-NY-NU
               endif
               if( c_str.ne.' ' )
     1         write(*,7010) k,name,index,'P',j, exact, approx,
     1              reldiff, c_str
            enddo
C            write(*,*)
            P(j) = P0(j)
         enddo

 1530    LAM(k) = 0.d0
      enddo

C
C     END
C
      stop


 7000 format(a3,'(',i4,') Eqn',i4,a3,3d17.8,a4)
 7010 format(i4,': row',a3,'(',i4,') col',a3,'(',i4,')',3d17.8,a4)
 8000 format(a)
 8010 format(i10)
 8020 format(3d23.15)

 9000 continue
      write(*,*) 'some maximal dimension exceeded - abort.'
      stop
C
C     That's it?
C
      end

      subroutine PERTURB(N, V, V0, VB, HV, name, EPS)
      implicit none
      integer N
      double precision V(N), V0(N), VB(2,N), HV(N), EPS
      character *(*) name
      integer i
      double precision eps2

      eps2 = sqrt(EPS)
C
C     Obtain slightly perturbed starting point (to get rid of inital zeros)
C
      do i = 1, N
         HV(i) = dmin1( 1.d0, VB(2,i)-VB(1,i) )
         if( HV(i).le.0.d0 ) then
            write(*,8100) name,i,VB(1,i),VB(2,i)
            stop
         endif
         if( V(i)+eps2*HV(i).ge.VB(2,i) ) then
            HV(i) = -HV(i)
         endif
         V0(i) = V(i) + eps2*HV(i)/2.d0
         HV(i) = HV(i)*EPS
      enddo

 8100 format('Error with bound for ',a,'(',i4,'): lo=',d23.15,
     1    ' up=',d23.15)

      return
      end
