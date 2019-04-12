C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: daemodel_f.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine DAEMODEL_F( NZ_OPT, NY_OPT, NU_OPT, NP_OPT, T, Z,
     1     ZDOT, Y, U_OPT, P_OPT, F )
C !DEC$ ATTRIBUTES DLLEXPORT :: DAEMODEL_F
C
C     Layer between decomposition code (dae2nlp) and actual model that
C     allows to give profiles to certain control variables and fixed values
C     to some parameters (see DYNAUX.INC)
C
C     Obtain residual of DAE model from MODEL_F
C
C     Author:   Andreas Waechter    09-05-01
C
      implicit none
C
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
      double precision F(NZ_OPT+NY_OPT) ! residual of DAE equations
C
      include 'DYNAUX.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/
C
      double precision U(NUMAX) ! controls given to model including given
                                ! profile
      double precision P(NPMAX) ! parameters given to model including fixed
                                ! ones

      integer l_uopt, l_popt, l_ele, i, GET_IELE4T
      double precision ratio, dt, dti

      NO_MODEL_F_CALL = NO_MODEL_F_CALL + 1

C     Only need to do something special if there are indeed given profiles 
C     or fixed parameters

      if( NU_PROF.eq.0 .and. NP_FIX.eq.0 ) then
         call MODEL_F( NZ, NY, NU, NP, T, Z, ZDOT, Y, U_OPT, P_OPT, F )
         return
      endif

C     Find element l_ele which is responsible for T
      l_ele = GET_IELE4T(T)
      if( l_ele.eq.NELE+1 ) l_ele = NELE

      dt  = T-TI(l_ele)
      dti = TI(l_ele+1)-TI(l_ele)

C     Copy relevant elements from U_OPT into U and include values for
C     given profiles (linear interpolations)
      l_uopt = 0
      do i = 1, NU
         if( IU_PROF(i).eq.0 ) then
            l_uopt = l_uopt + 1
            U(i) = U_OPT(l_uopt)
         else
            ratio = (U_PROF(i,l_ele+1) - U_PROF(i,l_ele))/dti
            U(i)  = U_PROF(i,l_ele) + ratio*dt
         endif
      enddo

C     Copy relevant elements from P_OPT into P and include values for
C     fixed parameters
      l_popt = 0
      do i = 1, NP
         if( IP_FIX(i).eq.0 ) then
            l_popt = l_popt + 1
            P(i) = P_OPT(l_popt)
         else
            P(i) = P_FIX(i)
         endif
      enddo

C     Call the model with augemented U and P
      call MODEL_F( NZ, NY, NU, NP, T, Z, ZDOT, Y, U, P, F )

C     That's it

      return
      end

