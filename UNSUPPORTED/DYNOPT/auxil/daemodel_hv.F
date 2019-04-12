C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: daemodel_hv.F 531 2004-03-11 01:31:07Z andreasw $
#ifdef OLDADOLC
      subroutine DAEMODEL_HV_OLD(NZ_OPT,NY_OPT,NU_OPT,NP_OPT,
     1     t,z,dz,y,u,p,LAM,V_Z,V_DZ,V_Y,V_U,V_P,
     1     HV_Z,HV_DZ,HV_Y,HV_U,HV_P)
      implicit none
      integer NZ_OPT, NY_OPT, NU_OPT, NP_OPT
      double precision t, z(nz_opt),dz(nz_opt),y(ny_opt)
      double precision u(nu_opt),p(np_opt)
      double precision LAM(nz_opt+ny_opt)
      double precision V_Z(nz_opt)
      double precision V_DZ(nz_opt)
      double precision V_Y(ny_opt)
      double precision V_U(nu_opt)
      double precision V_P(np_opt)
      double precision HV_Z(nz_opt)
      double precision HV_DZ(nz_opt)
      double precision HV_Y(ny_opt)
      double precision HV_U(nu_opt)
      double precision HV_P(np_opt)
      include 'DYNAUX.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/
      if( NU_PROF.eq.0 .and. NP_FIX.eq.0 ) then
         call MODEL_HV_OLD(NZ,NY,NU,NP,t,z,dz,y,u,p,LAM,
     1     V_Z,V_DZ,V_Y,V_U,V_P,
     1     HV_Z,HV_DZ,HV_Y,HV_U,HV_P)
         return
      else
         write(*,*) 'With old version ADOL-C only NU_PROF = 0 '//
     1'and NP_FIX = 0 implemented.'
         stop
      endif
      return
      end

#else
      subroutine DAEMODEL_HV( NRHS, NZ_OPT, NY_OPT, NU_OPT, NP_OPT,
     1     T, Z, ZDOT, Y, U_OPT, P_OPT, LAM, RHS_OPT, RES_OPT,
     1     LDRS_OPT)
C
C     Layer between decomposition code (dae2nlp) and actual model that
C     allows to give profiles to certain control variables and fixed values
C     to some parameters (see DYNAUX.INC)
C
C     Obtain Hessian-vector products from MODEL_HV
C
C     Author:   Andreas Waechter    09-05-01
C
      implicit none
C
      integer NRHS              ! number of right hand sides
      integer NZ_OPT            ! number of differential variables seen by
                                ! optimizer (is same as NZ in DYNAUX.INC)
      integer NY_OPT            ! number of algebraic variables seen by
                                ! optimizer (is same as NZ in DYNAUX.INC)
      integer NU_OPT            ! number of control variables seen by
                                ! optimizer
      integer NP_OPT            ! number of parameters seen by optimizer
      double precision T        ! time at which model is to be evaluated
      double precision Z(NZ_OPT) ! values of differenial variables
      double precision ZDOT(NZ_OPT) ! derivatives of differenial variables
      double precision Y(NY_OPT) ! values of algebraic variables
      double precision U_OPT(NU_OPT) ! values of control variables
      double precision P_OPT(NP_OPT) ! values of parameters
      double precision LAM(NZ_OPT+NY_OPT) ! multipliers for DAE equations
      integer LDRS_OPT          ! leading dimension of RHS_OPT and RES_OPT
      double precision RHS_OPT(LDRS_OPT, NRHS) ! vectors to be multiplied
      double precision RES_OPT(LDRS_OPT, NRHS) ! results
C
      include 'DYNAUX.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/
C
      double precision U(NUMAX) ! controls given to model including given
                                ! profile
      double precision P(NPMAX) ! parameters given to model including fixed
                                ! ones
      integer LDRS
      parameter( LDRS = 2*NZMAX+NYMAX+NUMAX+NPMAX )
      double precision RHS(LDRS, NCOLMAX+1)
                                ! right hand sides how model sees them
                                ! NOTE: cannot work with more than NCOLMAX
                                ! right hand sides
      double precision RES(LDRS, NCOLMAX+1)
                                ! results how model sees them
                                ! NOTE: cannot work with more than NCOLMAX
                                ! right hand sides

      integer l_uopt, l_popt, l_ele, l_resopt, i, GET_IELE4T
      double precision ratio, dt, dti

C     Only need to do something special if there are indeed given profiles 
C     or fixed parameters

      if( NU_PROF.eq.0 .and. NP_FIX.eq.0 ) then
         call MODEL_HV( NRHS, NZ, NY, NU, NP, T, Z, ZDOT, Y,
     1        U_OPT, P_OPT, LAM, RHS_OPT, RES_OPT, LDRS_OPT)
         return
      endif

C     Sanity check
      if( NRHS.gt.LDRS ) then
         call C_OUT(2,0,1,
     1         'Error in daemodel_hv: NRHS too large. Abort.')
         stop
      endif

C     Find element l_ele which is responsible for T
      l_ele = GET_IELE4T(T)
      if( l_ele.eq.NELE+1 ) l_ele = NELE

      do i = 1, 2*NZ+NY
         call DCOPY(NRHS, RES_OPT(i,1), LDRS_OPT, RES(i,1), LDRS)
      enddo

      dt  = T-TI(l_ele)
      dti = TI(l_ele+1)-TI(l_ele)

C     Copy relevant elements from U_OPT into U and include values for
C     given profiles (linear interpolations)
      l_uopt = 0
      do i = 1, NU
         if( IU_PROF(i).eq.0 ) then
            l_uopt = l_uopt + 1
            U(i) = U_OPT(l_uopt)
            call DCOPY(NRHS, RES_OPT(2*NZ+NY+l_uopt,1), LDRS_OPT,
     1           RES(2*NZ+NY+i,1), LDRS)
         else
            ratio = (U_PROF(i,l_ele+1) - U_PROF(i,l_ele))/dti
            U(i)  = U_PROF(i,l_ele) + ratio*dt
            call DCOPY(NRHS, 0.d0, 0, RES(2*NZ+NY+i,1), LDRS)
         endif
      enddo

C     Copy relevant elements from P_OPT into P and include values for
C     fixed parameters
      l_popt = 0
      do i = 1, NP
         if( IP_FIX(i).eq.0 ) then
            l_popt = l_popt + 1
            P(i) = P_OPT(l_popt)
            call DCOPY(NRHS, RES_OPT(2*NZ+NY+NU_OPT+l_popt,1), LDRS_OPT,
     1           RES(2*NZ+NY+NU+i,1), LDRS)
         else
            P(i) = P_FIX(i)
            call DCOPY(NRHS, 0.d0, 0, RES(2*NZ+NY+NU+i,1), LDRS)
         endif
      enddo

C     Call the model with augemented U and P
      call MODEL_HV( NRHS, NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P,
     1     LAM, RHS, RES, LDRS)

C     Copy RES from model into RES for optimizer
      do i = 1, 2*NZ+NY
         call DCOPY(NRHS, RES(i,1), LDRS,  RES_OPT(i,1), LDRS_OPT)
      enddo

      l_resopt = 2*NZ+NY
      do i = 1, NU
         if( IU_PROF(i).eq.0 ) then
            l_resopt = l_resopt + 1
            call DCOPY(NRHS, RES(2*NZ+NY+i,1), LDRS,
     1           RES_OPT(l_resopt,1), LDRS_OPT)
         endif
      enddo

      do i = 1, NP
         if( IP_FIX(i).eq.0 ) then
            l_resopt = l_resopt + 1
            call DCOPY(NRHS, RES(2*NZ+NY+NU+i,1), LDRS,
     1           RES_OPT(l_resopt,1), LDRS_OPT)
         endif
      enddo

C     That's it

      return
      end
#endif
