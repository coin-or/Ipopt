C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C     $Id: fprogs.F 531 2004-03-11 01:31:07Z andreasw $
c................................................................
      subroutine MODEL_F(NZ,NY,NU,NP,t,z,dz,y,u,p,F)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_F
      implicit none

      integer NZ, NY, NU, NP
      double precision t, z(NZ),dz(NZ),y(NY),u(NU),p(NP)
      double precision f(NZ+NY)

C     This is not standard f77:
      double precision xtmp(1+2*NZ+NY+NU+NP)
      integer retval, i, j, adolc_tag

#ifdef __MINGW32__
      integer FUNCTION
#elif WIN32
      interface
         integer function FUNCTION(TAG, M, N, ARG, RES)
C     !DEC$ ATTRIBUTES C, ALIAS:'_function_' :: FUNCTION
C     !DEC$ ATTRIBUTES REFERENCE:: TAG, M, N, ARG, RES
         integer TAG, M, N
         double precision ARG(*), RES(*)
         end function
      end interface
#else
      integer FUNCTION
#endif

      j = 1
      xtmp(j) = t

      if( nz.gt.0 ) then
         call DCOPY(nz, z, 1, xtmp(j+1), 1)
         j = j + nz
         call DCOPY(nz, dz, 1, xtmp(j+1), 1)
         j = j + nz
      endif

      if( ny.gt.0 ) then
         call DCOPY(ny, y, 1, xtmp(j+1), 1)
         j = j + ny
      endif

      if( nu.gt.0 ) then
         call DCOPY(nu, u, 1, xtmp(j+1), 1)
         j = j + nu
      endif

      if( np.gt.0 ) then
         call DCOPY(np, p, 1, xtmp(j+1), 1)
         j = j + np
      endif

      adolc_tag = 1
      retval = FUNCTION(adolc_tag, NZ+NY, j, xtmp, F)
      if( retval.lt.0 ) then
         write(*,*)
     1        'Error in MODEL_F: FUNCTION returns retval =', retval
         stop
      endif

      return
      end

c................................................................
      subroutine MODEL_DF(NZ,NY,NU,NP,t,z,dz,y,u,p,DF,LDDF)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_DF

      implicit none	

      integer NZ, NY, NU, NP, LDDF
      double precision t, z(nz),dz(nz),y(ny),u(nu),p(np)
      double precision df(LDDF,*)

C     This is not standard f77:
      double precision xtmp(1+2*nz+ny+nu+np)
      integer retval, j, adolc_tag
#ifdef __MINGW32__
      integer MYJAC
#elif WIN32
      interface
         integer function MYJAC(TAG, M, N, ARG, JACDIM, JAC)
C     !DEC$ ATTRIBUTES C, ALIAS:'_myjac_' :: MYJAC
C     !DEC$ ATTRIBUTES REFERENCE:: TAG, M, N, ARG, JACDIM, JAC
         integer TAG, M, N, JACDIM
         double precision ARG(*), JAC(*)
         end function
      end interface
#else
      integer MYJAC
#endif

      j = 1
      xtmp(j) = t

      if( nz.gt.0 ) then
         call DCOPY(nz, z, 1, xtmp(j+1), 1)
         j = j + nz
         call DCOPY(nz, dz, 1, xtmp(j+1), 1)
         j = j + nz
      endif

      if( ny.gt.0 ) then
         call DCOPY(ny, y, 1, xtmp(j+1), 1)
         j = j + ny
      endif

      if( nu.gt.0 ) then
         call DCOPY(nu, u, 1, xtmp(j+1), 1)
         j = j + nu
      endif

      if( np.gt.0 ) then
         call DCOPY(np, p, 1, xtmp(j+1), 1)
         j = j + np
      endif

C Note: Since we don't want dervatives with respect to t, in MYJAC
C       the first column is deleted!
      adolc_tag = 1
      retval = MYJAC(adolc_tag, nz+ny, j, xtmp, LDDF, df)
      if( retval.lt.0 ) then
         write(*,*) 'Error in MODEL_DF: MYJAC returns retval =', retval
         stop
      endif

      return
      end

c................................................................
      subroutine MODEL_DF_STRUC(NZ,NY,NU,NP,DF,LDDF)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_DF_STRUC

      implicit none
	
      integer NZ, NY, NU, NP, LDDF
      double precision df(LDDF,2*NZ+NY+NU+NP)

      integer retval, i, nnz, adolc_tag
#ifdef __MINGW32__
      integer MYJACPAT
#elif WIN32
      interface
         integer function MYJACPAT(TAG, M, N, NZ, AVAR, ACON)
C     !DEC$ ATTRIBUTES C, ALIAS:'_myjacpat_' :: MYJACPAT
C     !DEC$ ATTRIBUTES REFERENCE:: TAG, M, N, NZ, AVAR, ACON
         integer TAG, M, N, NZ, AVAR(*), ACON(*)
         end function
      end interface
#else 
      integer MYJACPAT
#endif
C
C     We might have to change the fill-in here in the
C
      integer NFILLIN
      parameter( NFILLIN = 20 )
      integer acon((NZ+NY)*NFILLIN), avar((NZ+NY)*NFILLIN)
C
C     Obtain sparcity pattern from ADOLC (in sparse MA28 notation)
C     (i.e. since this is what we want, we don't even need to store
C      everything in the big DF matrix...)
C
      nnz=(NZ+NY)*NFILLIN
      adolc_tag = 1
      retval = MYJACPAT(adolc_tag, nz+ny, 1+2*nz+ny+nu+np,
     1     nnz, avar, acon)
      if( retval.lt.0 ) then
         write(*,*)
     1        'Error in MODEL_DF_STRUC: MYJACPAT returns retval =',
     2        retval
         stop
      endif
C
C     Now fill DF
C
      do i = 1, 2*nz+ny+nu+np
         call DCOPY(nz+ny, 0d0, 0, DF(1,i), 1)
      enddo
C
C Note:  We don't want derivatives with respect to t
C
      do i = 1, nnz
         if( avar(i).gt.1 ) then
            DF(acon(i),avar(i)-1) = 1d0
         endif
      enddo

      return
      end

c................................................................
#ifdef OLDADOLC
      subroutine MODEL_HV_OLD(NZ,NY,NU,NP,t,z,dz,y,u,p,LAM,
     1     V_Z,V_DZ,V_Y,V_U,V_P,
     1     HV_Z,HV_DZ,HV_Y,HV_U,HV_P)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_HV_OLD
C
C     I just leave this in for possible debugging
C
      implicit none
	
      integer NZ, NY, NU, NP
      double precision t, z(nz),dz(nz),y(ny),u(nu),p(np)
      double precision LAM(nz+ny)
      double precision V_Z(nz)
      double precision V_DZ(nz)
      double precision V_Y(ny)
      double precision V_U(nu)
      double precision V_P(np)
      double precision HV_Z(nz)
      double precision HV_DZ(nz)
      double precision HV_Y(ny)
      double precision HV_U(nu)
      double precision HV_P(np)

C     This is not standard f77:
      double precision xtmp(1+2*nz+ny+nu+np)
      double precision vtmp(1+2*nz+ny+nu+np)
      double precision hvtmp(1+2*nz+ny+nu+np)
      integer retval, j, adolc_tag, i
#ifdef __MINGW32__
      integer LAGRA_HESS_VEC
#elif WIN32
      interface
         integer function LAGRA_HESS_VEC(TAG, M, N, X, V, LAM, HV)
C     !DEC$ ATTRIBUTES C, ALIAS:'_lagra_hess_vec_' :: LAGRA_HESS_VEC
C     !DEC$ ATTRIBUTES REFERENCE:: TAG, M, N, X, V, LAM, HV
         integer TAG, M, N
         double precision X(*), V(*), LAM(*), HV(*)
         end function
      end interface
#else
      integer LAGRA_HESS_VEC
#endif

      j = 1
      xtmp(j) = t
      vtmp(j) = 0.d0

      if( nz.gt.0 ) then
         call DCOPY(nz,   z, 1, xtmp(j+1), 1)
         call DCOPY(nz, v_z, 1, vtmp(j+1), 1)
         j = j + nz
         call DCOPY(nz,   dz, 1, xtmp(j+1), 1)
         call DCOPY(nz, v_dz, 1, vtmp(j+1), 1)
         j = j + nz
      endif

      if( ny.gt.0 ) then
         call DCOPY(ny,   y, 1, xtmp(j+1), 1)
         call DCOPY(ny, v_y, 1, vtmp(j+1), 1)
         j = j + ny
      endif

      if( nu.gt.0 ) then
         call DCOPY(nu,   u, 1, xtmp(j+1), 1)
         call DCOPY(nu, v_u, 1, vtmp(j+1), 1)
         j = j + nu
      endif

      if( np.gt.0 ) then
         call DCOPY(np,   p, 1, xtmp(j+1), 1)
         call DCOPY(np, v_p, 1, vtmp(j+1), 1)
         j = j + np
      endif

      adolc_tag = 1
      retval = LAGRA_HESS_VEC(adolc_tag, nz+ny, j, xtmp, vtmp, LAM,
     1     hvtmp)
      if( retval.lt.0 ) then
         write(*,*)
     1        'Error in MODEL_DF: LAGRA_HESS_VEC returns retval =',
     2        retval
         stop
      endif

      j = 1
      if( nz.gt.0 ) then
         call DCOPY(nz, hvtmp(j+1), 1, hv_z, 1)
         j = j + nz
         call DCOPY(nz, hvtmp(j+1), 1, hv_dz, 1)
         j = j + nz
      endif

      if( ny.gt.0 ) then
         call DCOPY(ny, hvtmp(j+1), 1, hv_y, 1)
         j = j + ny
      endif

      if( nu.gt.0 ) then
         call DCOPY(nu, hvtmp(j+1), 1, hv_u, 1)
         j = j + nu
      endif

      if( np.gt.0 ) then
         call DCOPY(np, hvtmp(j+1), 1, hv_p, 1)
         j = j + np
      endif
      
      return
      end
c................................................................
#else
      subroutine MODEL_HV(NRHS,NZ,NY,NU,NP,t,z,dz,y,u,p,
     1     LAM,RHS,RES,LDRS)
C !DEC$ ATTRIBUTES DLLEXPORT :: MODEL_HV

      implicit none
	
      integer NRHS, NZ, NY, NU, NP
      double precision t, z(nz),dz(nz),y(ny),u(nu),p(np)
      double precision LAM(nz+ny)
      integer LDRS
      double precision RHS(LDRS, NRHS)
      double precision RES(LDRS, NRHS)

C     This is not standard f77:
      double precision xtmp(1+2*nz+ny+nu+np)

      integer retval, j, adolc_tag, i
#ifdef __MINGW32__
      integer MYHESSMAT
#elif WIN32
      interface
         integer function MYHESSMAT(TAG, M, N, NRHS, X, LAM,
     1        RHS, RES, LDRS)
C     !DEC$ ATTRIBUTES C, ALIAS:'_myhessmat_' :: MYHESSMAT
C     !DEC$ ATTRIBUTES REFERENCE:: TAG, M, N, NRHS, X, LAM
C     !DEC$ ATTRIBUTES REFERENCE:: RHS, RES, LDRS
         integer TAG, M, N, NRHS, LDRS
         double precision X(*), LAM(*), RHS(*), RES(*)
         end function
      end interface
#else
      integer MYHESSMAT
#endif

      j = 1
      xtmp(j) = t

      if( nz.gt.0 ) then
         call DCOPY(nz,   z, 1, xtmp(j+1), 1)
         j = j + nz
         call DCOPY(nz,  dz, 1, xtmp(j+1), 1)
         j = j + nz
      endif

      if( ny.gt.0 ) then
         call DCOPY(ny,   y, 1, xtmp(j+1), 1)
         j = j + ny
      endif

      if( nu.gt.0 ) then
         call DCOPY(nu,   u, 1, xtmp(j+1), 1)
         j = j + nu
      endif

      if( np.gt.0 ) then
         call DCOPY(np,   p, 1, xtmp(j+1), 1)
         j = j + np
      endif

      adolc_tag = 1
      retval = MYHESSMAT(adolc_tag, nz+ny, j, NRHS, xtmp, LAM,
     1     RHS, RES, LDRS)
      if( retval.lt.0 ) then
         write(*,*)
     1        'Error in MODEL_DF: MYHESSMAT returns retval =',
     2        retval
         stop
      endif

      return
      end
#endif
