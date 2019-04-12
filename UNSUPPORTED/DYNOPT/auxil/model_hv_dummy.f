C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: model_hv_dummy.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine MODEL_HV( NRHS, NZ_OPT, NY_OPT, NU_OPT, NP_OPT,
     1     T, Z, ZDOT, Y, U_OPT, P_OPT, LAM, RHS_OPT, RES_OPT,
     1     LDRS_OPT)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_HV
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

      write(*,*) 'Error: Shouldn''t get into this MODEL_HV!'
      stop

      return
      end

