C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: DYNOPC.INC 531 2004-03-11 01:31:07Z andreasw $
C
C     Parameters and common blocks for DYNOPC
C
C     Authors:   Yi-dong Lang
C                Andreas Waechter    09-03-01   cleaned up
C

C     Information for objective function
      integer NOBJMAX
      parameter( NOBJMAX = 10 )
      integer NZ_IN_OBJ, NY_IN_OBJ, NP_IN_OBJ, NE_IN_OBJ
                                ! number of variable in objective function
      integer Z_INDEX(NOBJMAX), Y_INDEX(NOBJMAX)
      integer E_INDEX(NOBJMAX), P_INDEX(NOBJMAX)
      double precision AIJZ(NOBJMAX,NOBJMAX)
      double precision BIJY(NOBJMAX,NOBJMAX), CIP(NOBJMAX)
C     IGNORE:
C      integer obj_mode          ! mode of objective function 
C                                !   = 0  calculated at final time
C                                !   = 1  calculated at specified elements


C     Unit number for C_OUT output in IPOPT
      integer OutputUnit_COUT
      
C     Final time
      double precision end_time

C----- for additional bound
      integer IADB              ! nonzero, if there are
                                ! additional bound constraints
      integer NADBMAX
      parameter( NADBMAX = 5 )
      integer nz_in_adb, ny_in_adb, nu_in_adb, ne_in_adb
      integer attrz(NADBMAX,NADBMAX)
      integer attry(NADBMAX,NADBMAX)
      integer attru(NADBMAX,NADBMAX) 
      integer indzadb(NADBMAX), indyadb(NADBMAX)
      integer induadb(NADBMAX), indeadb(NADBMAX)
      double precision z_adb(NADBMAX,NADBMAX)
      double precision y_adb(NADBMAX,NADBMAX)
      double precision u_adb(NADBMAX,NADBMAX)
      integer nzadb, nyadb, nuadb ! Counter for added bounds
      integer izadb(NADBMAX), iyadb(NADBMAX), iuadb(NADBMAX)
                                ! array to store index of i
      integer jzadb(NADBMAX), jyadb(NADBMAX), juadb(NADBMAX)
                                ! array to store index of j
      integer iele_adb(NADBMAX)  ! the element that timepoint is in
      double precision time_adb(NADBMAX)

C     Some graphics stuff...
      integer NDSPLMAX
      parameter( NDSPLMAX = 12 )
      integer indexz(NDSPLMAX), indexy(NDSPLMAX),indexu(NDSPLMAX) 
                                        ! the index of variables to be 
                                        ! displayed
      integer ndsplz, ndsply, ndsplu    ! numbers of variable that user want to
                                        ! get display

C     Flag for JAC_NLEQ1S - set to .FALSE. before first call of NLEQ1S
C     in initcond.f
      logical JAC_INIT

C     Want this out if not OLDSTYLE
      double precision eps3

      common /DYNOPC/ end_time, AIJZ, BIJY, CIP,
     1     z_adb, y_adb, u_adb, time_adb,
     1     eps3,                ! This is only for OLDSTYLE version
     1     NZ_IN_OBJ, NY_IN_OBJ, NP_IN_OBJ, NE_IN_OBJ,
     1     Z_INDEX, Y_INDEX, E_INDEX, P_INDEX,
     1     IADB, ndsplz, ndsply, ndsplu,
     1     indexz, indexy, indexu,
     1     OutputUnit_COUT, JAC_INIT,
     1     nz_in_adb, ny_in_adb, nu_in_adb, ne_in_adb,
     1     attrz, attry, attru, indzadb, indyadb, induadb, indeadb,
     1     nzadb, nyadb, nuadb,
     1     izadb, iyadb, iuadb, jzadb, jyadb, juadb,
     1     iele_adb
 
