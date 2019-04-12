C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine REORDER_IB(N, NORIG, NLB, ILB, NUB, IUB,
     1                      IVAROLD, IVARNEW, ivarnew1, ib)
C
C*******************************************************************************
C
C    $Id: reorder_ib.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Reorder ILB and IUB from partition in IVAROLD to IVARNEW
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
CA    Andreas Waechter      05/01/02  Release as version IPOPT 2.0
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
CP   N         I    INT    number of variables (without fixed vars)
CP   NORIG     I    INT    total number of variables (including fixed)
CP   NLB       I    INT    number of lower bounds (no fixing entries)
CP   ILB      I/O   INT    BNDS_L(i) corresponds to X(ILB(i))
CP                            I: ordered like IVAROLD
CP                            O: ordered like IVARNEW
CP   NUB       I    INT    number of upper bounds (no fixing entries)
CP   IUB      I/O   INT    BNDS_U(i) corresponds to X(IUB(i))
CP                            I: ordered like IVAROLD
CP                            O: ordered like IVARNEW
CP   IVAROLD   O    INT    information about old partitioning
CP                            i = 1..M      XORIG(IVAROLD(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAROLD(i)) independent
CP                            Note: fixed variables do not occur in IVAROLD
CP   IVARNEW   O    INT    information about old partitioning
CP                            i = 1..M      XORIG(IVARNEW(i)) dependent
CP                            i = (M+1)..N  XORIG(IVARNEW(i)) independent
CP                            Note: fixed variables do not occur in IVARNEW
CP   IVARNEW1  O    INT    inverse of IVARNEW
CP   ib        W
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
CS
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

C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer N
      integer NORIG
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      integer IVAROLD(N)
      integer IVARNEW(N)
      integer ivarnew1(NORIG)
      integer ib(N)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Compute inverse of IVARNEW
C
      if( NORIG.gt.N ) then
         do i = 1, NORIG
            ivarnew1(i) = 0
         enddo
      endif
      do i = 1, N
         ivarnew1(IVARNEW(i)) = i
      enddo
C
C     Change ILB
C
      do i = 1, NLB
         ib(i) = ivarnew1(IVAROLD(ILB(i)))
      enddo
      do i = 1, NLB
         ILB(i) = ib(i)
      enddo
C
C     Change IUB
C
      do i = 1, NUB
         ib(i) = ivarnew1(IVAROLD(IUB(i)))
      enddo
      do i = 1, NUB
         IUB(i) = ib(i)
      enddo

      return
      end
