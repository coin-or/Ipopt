C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

      subroutine draw_iter_info(iter,x,rg,infsbl,
     1                       alfaV,alfaTau,alfaX,kkt,obj,mu)

      USE dflib
         
      implicit none

      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/

      integer(4) i4
      logical(4) results
      logical ex
      integer iter
	double precision x(*)
      double precision rg, infsbl, alfaV, alfaTau, alfaX, kkt, obj, mu
      character*80 fname

      i4 =setactiveqq(CurveUnit)
      call draw_profile(1,x)
 
      call draw_opt_path(alfaV,alfaTau,alfaX, kkt, obj,iter)

      call draw_rg_inf_path(rg,infsbl,mu)

     
C      if (iter .eq. 1) then
C	   call GET_FILENAME(
C         fname = trim(fname_root0)
C         inquire(file=trim(fname), exist = ex)
C         if(.not. ex) then
C            results = Makedirqq(trim(fname))  
C         end if       
C          fname = trim(fname_root0)//"\history.dat"   
C         open(31,file = trim(fname),status='unknown')
C         write(31, 100)
C      end if
C      write(31,110) iter, kkt, obj, alfaV, alfaTau, alfaX, rg, infsbl

C      if(kkt .lt. opt_eps) close(31)      

100   format(1x,
     !"ITER       KKT       OBJ      AlfaV  AlfaTau  AlfaX
     !    RG     INFSBL "/)
110   format (1x, i3,4x, e9.3e1,2x,e9.3e1, 3(3x, f5.3), 2(2x e7.2e2))

      return    
      end

C
C     This is for ITER_OUT_DYNOPC.F
C
      subroutine SET_DRAWSUB(DRAWSUB)
C
C     Setting the common block variable
C
      implicit none
      integer(4) DRAWSUB
      include 'DYNAUX.INC'
	include 'DYNOPC.INC'
      include 'DYNGRA.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/, /GRAPH/
      DRAWSUBP = %LOC(DRAWSUB)
      return
      end
