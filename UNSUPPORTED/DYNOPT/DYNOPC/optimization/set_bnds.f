C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: set_bnds.f 531 2004-03-11 01:31:07Z andreasw $
      subroutine SET_BNDS(ZB, YB, UB, PB, NLB, NUB, ILB, IUB,
     1     BNDS_L, BNDS_U)
C
C     Set the bound information in IPOPT style
C
C     Authors:  Yidong Lang, Andreas Waechter   10-02-01
C
      implicit none
C
      double precision ZB(2,*), YB(2,*), UB(2,*), PB(2,*)
      integer NLB, NUB, ILB(*), IUB(*)
      double precision BNDS_L(*), BNDS_U(*)
C
      include 'DYNAUX.INC'
      include 'DYNOPC.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/, /DYNOPC/
C
      integer i, lv, lx, iele, icol
      integer llbz, lubz, llby, luby, llbu, lubu

C     First set the regular bounds in each element

      NLB = 0
      NUB = 0
C
C     Differential variables
C
      do i = 1, NZ
         if( ZB(1,i).gt.-VLARGE ) then
            lx = NZ + NCOL*(NZ+NY) + i ! No bounds on initial conditions
            do iele = 1, NELE
               NLB         = NLB + 1
               ILB(NLB)    = lx
               BNDS_L(NLB) = ZB(1,i)
               lx = lx + NZ + NCOL*(NZ+NY)
            enddo
         endif
         if( ZB(2,i).lt. VLARGE ) then
            lx = NZ + NCOL*(NZ+NY) + i ! No bounds on initial conditions
            do iele = 1, NELE
               NUB         = NUB + 1
               IUB(NUB)    = lx
               BNDS_U(NUB) = ZB(2,i)
               lx = lx + NZ + NCOL*(NZ+NY)
            enddo
         endif
      enddo
      llbz = NLB
      lubz = NUB
C
C     Algebraic variables
C
      do i = 1, NY
         if( YB(1,i).gt.-VLARGE ) then
            lx = 2*NZ + i
            do iele = 1, NELE
               do icol = 1, NCOL
                  NLB         = NLB + 1
                  ILB(NLB)    = lx
                  BNDS_L(NLB) = YB(1,i)
                  lx = lx + (NZ+NY)
               enddo
               lx = lx + NZ
            enddo
         endif
         if( YB(2,i).lt. VLARGE ) then
            lx = 2*NZ + i
            do iele = 1, NELE
               do icol = 1, NCOL
                  NUB         = NUB + 1
                  IUB(NUB)    = lx
                  BNDS_U(NUB) = YB(2,i)
                  lx = lx + (NZ+NY)
               enddo
               lx = lx + NZ
            enddo
         endif
      enddo
      llby = NLB
      luby = NUB
C
C     Control variables
C
      do i = 1, NU
         if( UB(1,i).gt.-VLARGE ) then
            lx = NZ + NELE*(NCOL*(NZ+NY)+NZ) + i
            do iele = 1, NELE
               do icol = 1, NCOL
                  NLB         = NLB + 1
                  ILB(NLB)    = lx
                  BNDS_L(NLB) = UB(1,i)
                  lx = lx + NU
               enddo
            enddo
         endif
         if( UB(2,i).lt. VLARGE ) then
            lx = NZ + NELE*(NCOL*(NZ+NY)+NZ) + i
            do iele = 1, NELE
               do icol = 1, NCOL
                  NUB         = NUB + 1
                  IUB(NUB)    = lx
                  BNDS_U(NUB) = UB(2,i)
                  lx = lx + NU
               enddo
            enddo
         endif
      enddo
      llbu = NLB
      lubu = NUB
C
C     Parameters
C
      lx = NZ + NELE*(NCOL*(NZ+NY+NU)+NZ)
      do i = 1, NP
         if( PB(1,i).gt.-VLARGE ) then
            NLB         = NLB + 1
            ILB(NLB)    = lx + i
            BNDS_L(NLB) = PB(1,i)
         endif
         if( PB(2,i).lt. VLARGE ) then
            NUB         = NUB + 1
            IUB(NUB)    = lx + i
            BNDS_U(NUB) = PB(2,i)
         endif
      enddo
C
C     Now let's take care of additional bounds
C
      if( IADB.ne.0 ) then

C     Differential variables

         do i = 1, NZADB
            lv = INDZADB(IZADB(i))
            lx = IELE_ADB(JZADB(i))*(NCOL*(NZ+NY)+NZ) + lv
            call SET_ADB(ATTRZ(IZADB(i),JZADB(i)), lv, lx,
     1           Z_ADB(IZADB(i),JZADB(i)), ZB, NLB, NUB, 1, llbz,
     1           1, lubz, ILB, IUB, BNDS_L, BNDS_U)
         enddo

C     Algebraic variables

         do i = 1, NYADB
            lv = INDYADB(IYADB(i))
            lx = 2*NZ + (IELE_ADB(JYADB(i))-1)*(NCOL*(NZ+NY)+NZ) + lv
c            do icol = 0, NCOL-1
            do icol = NCOL-1, NCOL-1
               call SET_ADB(ATTRY(IYADB(i),JYADB(i)), lv,
     1              lx+icol*(NZ+NY), Y_ADB(IYADB(i),JYADB(i)), YB,
     1              NLB, NUB, llbz+1, llby, lubz+1, luby, ILB, IUB,
     1              BNDS_L, BNDS_U)
            enddo
         enddo

C     Control variables

         do i = 1, NUADB
            lv = INDUADB(IUADB(i))
            lx = NZ + NELE*(NCOL*(NZ+NY)+NZ) +
     1           (IELE_ADB(JYADB(i))-1)*NCOL*NU + lv
            do icol = 0, NCOL-1
               call SET_ADB(ATTRU(IUADB(i),JUADB(i)), lv,
     1              lx+icol*NU, U_ADB(IUADB(i),JUADB(i)), UB,
     1              NLB, NUB, llby+1, llbu, luby+1, lubu, ILB, IUB,
     1              BNDS_L, BNDS_U)
            enddo
         enddo
      endif

C     That's it!

C
CDELETEME
      if( .false. ) then
         write(*,*) 'BOUNDS FOR TE included!'
         do i = 1, NZ
            nlb = nlb + 1
            bnds_l(nlb) = - 1e-6 *max(1,abs(ZINIT(i)))
            ilb(nlb) = NZ + (NELE-1)*(NCOL*(NZ+NY)+NZ)
     1           + (NCOL-1)*(NZ+NY) + i
            nub = nub + 1
            bnds_u(nub) = + 1e-6 *max(1,abs(ZINIT(i)))
            iub(nub) = ilb(nlb)
         enddo
      endif





      return
      end

C
C     auxilliary subroutine for assigning additional bounds
C
      subroutine SET_ADB(ATTR, LV, LX, V_ADB, VB, NLB, NUB,
     1     LLB1, LLB2, LUB1, LUB2, ILB, IUB, BNDS_L, BNDS_U)
C
      implicit none
      integer ATTR, LV, LX
      double precision V_ADB, VB(2,*)
      integer NLB, NUB, LLB1, LLB2, LUB1, LUB2
      integer ILB(*), IUB(*)
      double precision BNDS_L(*), BNDS_U(*)

      include 'DYNAUX.INC'
!DEC$ ATTRIBUTES DLLIMPORT :: /DYNAUX/

      integer i

      if( ATTR.eq.0 ) then      ! lower bound
         if( VB(1,LV).gt.-VLARGE ) then ! has been set before
            do i = LLB1, LLB2
               if( ILB(i).eq.LX ) then
                  BNDS_L(i) = V_ADB
                  return
               endif
            enddo
            write(*,*) 'ERROR in SET_ADB!' ! Take out of code later
            stop
         else
            NLB         = NLB + 1
            ILB(NLB)    = LX
            BNDS_L(NLB) = V_ADB
         endif
      else                      ! upper bound
         if( VB(2,LV).lt. VLARGE ) then ! has been set before
            do i = LUB1, LUB2
               if( IUB(i).eq.LX ) then
                  BNDS_U(i) = V_ADB
                  return
               endif
            enddo
            write(*,*) 'ERROR in SET_ADB!' ! Take out of code later
            stop
         else
            NUB         = NUB + 1
            IUB(NUB)    = LX
            BNDS_U(NUB) = V_ADB
         endif
      endif
      return
      end
