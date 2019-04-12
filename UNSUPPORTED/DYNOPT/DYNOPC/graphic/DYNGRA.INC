C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C     I guess, these are for Windows

!/ index of eidt box in initial group
      integer, parameter :: idb = 1026
      integer, parameter :: idb_cap = 1043 
      integer, parameter :: idb_b = 1063  ! index of edit_b1 in bounds
                                          ! dialog box
      integer, parameter :: id_vid = 1150 ! index of IDC_index1 in objfunc
      integer, parameter :: id_vc  = 1160 ! index of IDC_coef1 in objfunc
      integer, parameter :: id_pid = 1174
!/----------------------------------------------------
!      index of the edit boxes
!/----------------------------------------------------
      integer, parameter:: id_z  = 1085, id_a =  1095, id_ik = 1121,
     &                     id_y  = 1105, id_b =  1110, id_jk = 1131,
     &                     id_u  = 1136, id_c =  1138,
     &                     id_pi = 1140, id_pj = 1144

      logical   retlog
      logical*4 retlog4
      integer   retint
      integer*4 retint4, i_25

!/ counter of next button press
      integer   inext

      integer unit_w, unit_h     ! the unit size of a character in term of pixel.
      double precision ratiox, ratioy !ratio of in unit of pixels to character
      double precision origin, spanx, spany, factorx, factory 
      double precision final_value
      double precision xx0,yy0,xx1,yy1
      double precision temp0, opt_eps

      logical opt_flag ! set to true if this is the graphic for optimization

      integer maxx, maxy

C      integer,parameter:: 
C     1                       OutputUnit_1 = 10 ,! for coldae output
C     2                       OutputUnit_2 = 11 ,! for OPT output
C     3                       CurveUnit    = 12 ,! for variable iteration profile
C     4                       OutputUnit_0 = 9   ! for DDASSL (and all IPOPT C_OUT output)

      integer,parameter::
     4                       OutputUnit_0 = 21 ,! for DDASSL (and all IPOPT C_OUT output)
     1                       OutputUnit_1 = 10 ,! for coldae output
     2                       OutputUnit_2 = 11 ,! for OPT output
     3                       CurveUnit    = 12 ,! for OPT iteration profiles
     4                       CoverUnit    = 13 ,! for cover page
     5                       ObjecF_W     = 14 ,! for obj values vs ITER No.
     6                       YpyZpz_W     = 15 ,! for OPT path
     7                       Report_W     = 16 ,! for final OPT solution profiles
     8                       opt_path_w   = 17 ,
     9                       ypy_zpz_w    = 18
      integer,parameter::
     1                       fispace_Unit = 19 ,     !for read and write fispace.d
     2                       output_unit = 20

      integer,parameter::
     1                       nuinit_Unit  = 1  ,! for reading data from nuinit.d file
     2                       nucoef_Unit = 2    ,! for reading data from nucoef.d file
     3                       lastx_Unit   = 3   ,! for reading data from latex.d file
     4                       ct210_Unit   = 7  ,! for reading data from ct210.d file
     5                       nucomp_Unit = 8    ,
     6                       nucomp1_unit      = 9
       integer,parameter::
     1                       start_dat_U  = 4   ,! for writing data to startp.d file
     2                       output_dat_U = 5    ! for writing data to output.d file


      integer PRINTSUBP, DRAWSUBP    ! very dirty programming :)


C  The following are lower and upper bounds for the plotted variables
      double precision ZBG(2,NZMAX)
      double precision YBG(2,NYMAX)
      double precision UBG(2,NUMAX)
      double precision PBG(2,NPMAX)

      common /GRAPH/ ZBG, YBG, UBG, PBG,
     1  origin, spanx, spany, factorx, factory, final_value,
     1  xx0,yy0,xx1,yy1, temp0, ratiox, ratioy, opt_eps,
     1  maxx, maxy, unit_w, unit_h, opt_flag,
     1  retlog, retlog4, retint, retint4, i_25, inext, PRINTSUBP,
	1  DRAWSUBP
      save /GRAPH/
