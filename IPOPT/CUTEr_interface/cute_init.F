C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine CUTE_INIT(N, M, LX, X, NLB, ILB, BNDS_L,
     1                     NUB, IUB, BNDS_U)
C
C*******************************************************************************
C
C    $Id: cute_init.F 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Initialize interface to CUTE problem
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB
C
C-------------------------------------------------------------------------------
C                             Author, date
C-------------------------------------------------------------------------------
C
CA    Andreas Waechter      02/25/99
CA    Andreas Waechter      07/19/99  initialize slacks based on C
C
C-------------------------------------------------------------------------------
C                             Documentation
C-------------------------------------------------------------------------------
C
CD
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
C    Name     I/O   Type   Meaning
C
CP   N         O    INT    number of variables in problem statement
CP                            (including slacks for inequality constraints)
CP   M         O    INT    number of equality constraints
CP   LX        I    INT    actual declared length of X, LIB, BNDS_L, IUB, BNDS_U
CP   X         O    DP     starting point
CP   NLB       O    INT    number of lower bounds
CP   ILB       O    INT    indices for lower bounds
CP                            ( BNDS_L(i) is lower bound for X(ILB(i)) )
CP   BNDS_L    O    INT    values of lower bounds
CP   NUB       O    INT    number of upper bounds
CP   IUB       O    INT    indices for upper bounds
CP                            ( BNDS_U(i) is lower bound for X(IUB(i)) )
CP   BNDS_U    O    INT    values of upper bounds
C
C-------------------------------------------------------------------------------
C                             local variables
C-------------------------------------------------------------------------------
C
CL
C
C-------------------------------------------------------------------------------
C                             used subroutines
C-------------------------------------------------------------------------------
C
CS    CSTEUP
C
C*******************************************************************************
C
C                              Declarations
C
C*******************************************************************************
C
      IMPLICIT NONE
C
C*******************************************************************************
C
C                              Include files
C
C*******************************************************************************
C
#include "CUTE.INC"
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer M
      integer LX
      double precision X(LX)
      integer NLB
      integer ILB(LX)
      double precision BNDS_L(LX)
      integer NUB
      integer IUB(LX)
      double precision BNDS_U(LX)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision bl(CUTE_NMAX), bu(CUTE_NMAX)
      double precision v(CUTE_MMAX), cl(CUTE_MMAX), cu(CUTE_MMAX)
      double precision c(CUTE_MMAX), f, dummy
      logical equatn(CUTE_MMAX), linear(CUTE_MMAX)
      integer cnr_input, iout, i
      logical efirst, lfirst, nvfrst, ex
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
C-------------------------------------------------------------------------------

C
C     Call CSETUP to obtain problem size and starting point
C
      cnr_input = 60
      iout = 6
      efirst = .false.
      lfirst = .false.
      nvfrst = .false.

      open(cnr_input,file='OUTSDIF.d',status='old')

      call CSETUP(cnr_input, iout, CUTE_N, CUTE_M,
     1     X, bl, bu, CUTE_NMAX, equatn, linear, v, cl, cu,
     2     CUTE_MMAX, efirst, lfirst, nvfrst)
      close(cnr_input)
C
C     Added this:  Compute C in order to initialize slacks better
C
      call CCFG(CUTE_N, CUTE_M, X, M, C, .false., 1, 1, dummy, .false.)
C
      M = CUTE_M
      N = CUTE_N
C
C     Obtain bounds on variables
C
      NLB = 0
      do i = 1, CUTE_N
         if( bl(i).gt.CUTE_NOLB ) then
            NLB = NLB + 1
            ILB(NLB) = i
            BNDS_L(NLB) = bl(i)
         endif
      enddo

      NUB = 0
      do i = 1, CUTE_N
         if( bu(i).lt.CUTE_NOUB ) then
            NUB = NUB + 1
            IUB(NUB) = i
            BNDS_U(NUB) = bu(i)
         endif
      enddo
C
C     Find inequalities and augment X
C
      CUTE_NIQ = 0
      CUTE_NEQ = 0
      do i = 1, CUTE_M
         if( .not.equatn(i) ) then
            CUTE_NIQ = CUTE_NIQ + 1
            CUTE_IIQ(CUTE_NIQ) = i
            N = N + 1
            if( N.gt.LX ) then
               write(6,*) 'Error in cute_init: LX = ',
     1                    LX,' is too small. Abort.'
               stop
            endif
C
C           This is kind of arbitrary: initialize slack to zero...
C
C            X(N) = 0.d0
            X(N) = C(i)
            if( cl(i).gt.CUTE_NOLB ) then
               NLB = NLB + 1
               ILB(NLB) = N
               BNDS_L(NLB) = cl(i)
            endif
            if( cu(i).lt.CUTE_NOUB ) then
               NUB = NUB + 1
               IUB(NUB) = N
               BNDS_U(NUB) = cu(i)
            endif
         else
            CUTE_NEQ = CUTE_NEQ + 1
            CUTE_CEQ(CUTE_NEQ) = cl(i)
         endif
      enddo
C
C     For basis selection, write indices of slack variables into file
C     SLACKS.DAT, or delete this file, if no slack variables
C
      if( CUTE_NIQ.gt.0 ) then
         open(10,file='SLACKS.DAT',status='unknown')
         do i = CUTE_N+1, CUTE_N+CUTE_NIQ
            write(10,1000) i
 1000       format(i10)
         enddo
         close(10)
      else
         inquire(file='SLACKS.DAT',exist=ex)
         if( ex ) then
            open(10,file='SLACKS.DAT',status='old')
            close(10,status='delete')
         endif
      endif

 9999 continue
      return
      end
