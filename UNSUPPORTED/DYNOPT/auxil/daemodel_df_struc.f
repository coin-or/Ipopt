C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: daemodel_df_struc.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine DAEMODEL_DF_STRUC( NZ_OPT, NY_OPT, NU_OPT, NP_OPT,
     1     DF_OPT, LDDF )
C !DEC$ ATTRIBUTES DLLEXPORT :: DAEMODEL_DF_STRUC
C
C     Layer between decomposition code (dae2nlp) and actual model that
C     allows to give profiles to certain control variables and fixed values
C     to some parameters (see DYNAUX.INC)
C
C     Obtain structure of Jacobian of DAE model from MODEL_DF_STRUC
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
      integer LDDF              ! leading dimension of DF
      double precision DF_OPT(LDDF,2*NZ_OPT+NY_OPT+NU_OPT+NP_OPT)
                                ! Jacobian of DAEs how optimizer needs it
C
      include 'DYNAUX.INC'
C !DEC$ ATTRIBUTES DLLEXPORT :: /DYNAUX/
C
      double precision DF(NZMAX+NYMAX,2*NZMAX+NYMAX+NUMAX+NPMAX)
                                ! Jacobian of DAEs how models sees it
                                ! NOTE: This is not very memory efficient...

      integer l_dfopt, i
      double precision ratio, dt, dti

C     Only need to do something special if there are indeed given profiles 
C     or fixed parameters

      if( NU_PROF.eq.0 .and. NP_FIX.eq.0 ) then
         call MODEL_DF_STRUC( NZ, NY, NU, NP, DF_OPT, LDDF )
         return
      endif

C     Call the model to get structure for larger Jacobian after initializing
C     DF to zero
      do i = 1, 2*NZ+NY+NU+NP
         call DCOPY(NZ+NY, 0.d0, 0, DF(1,i), 1)
      enddo
      call MODEL_DF_STRUC( NZ, NY, NU, NP, DF, NZMAX+NYMAX )

C     Copy Jacobian from model into Jacobian for optimizer
      do i = 1, 2*NZ+NY
         call DCOPY(NZ+NY, DF(1,i), 1, DF_OPT(1,i), 1)
      enddo

      l_dfopt = 2*NZ+NY
      do i = 1, NU
         if( IU_PROF(i).eq.0 ) then
            l_dfopt = l_dfopt + 1
            call DCOPY(NZ+NY, DF(1,2*NZ+NY+i), 1, DF_OPT(1, l_dfopt), 1)
         endif
      enddo

      do i = 1, NP
         if( IP_FIX(i).eq.0 ) then
            l_dfopt = l_dfopt + 1
            call DCOPY(NZ+NY, DF(1,2*NZ+NY+NU+i), 1,
     1           DF_OPT(1, l_dfopt), 1)
         endif
      enddo

C     That's it

      return
      end

