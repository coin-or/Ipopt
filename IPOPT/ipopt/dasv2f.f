C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine DASV2F(N, NLB, ILB, NUB, IUB, SP_L, SP_U, FULL)
C
C*******************************************************************************
C
C    $Id: dasv2f.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Add two sparse vectors to full
CT     The whole thing could be improved by making FULL also sparse!
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
CP   N         I    INT    number of variables
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP                            SP_L(i) corresponds to FULL(ILB(i))
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            SP_U(i) corresponds to FULL(LUB(i))
CP   SP_L      I    DP     sparse vector for lower bounds
CP   SP_U      I    DP     sparse vector for upper bounds
CP   FULL      O    DP     full vector with sum of both sparse vectors
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
CS    DCOPY
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
      integer NLB
      integer ILB(NLB)
      integer NUB
      integer IUB(NUB)
      double precision SP_L(NLB)
      double precision SP_U(NUB)
      double precision FULL(N)
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
      call DCOPY(N, 0d0, 0, FULL, 1)
      do i = 1, NLB
         FULL(ILB(i)) = SP_L(i)
      enddo
      do i = 1, NUB
         FULL(IUB(i)) = FULL(IUB(i)) + SP_U(i)
      enddo

      return
      end
