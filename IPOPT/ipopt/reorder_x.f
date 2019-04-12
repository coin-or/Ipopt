C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine REORDER_X(N, NORIG, X, IVAROLD, IVARNEW, xtmp)
C
C*******************************************************************************
C
C    $Id: reorder_x.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Reorder X (or G) from partition in IVAROLD to IVARNEW
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
CP   X        I/O   DP     I: partitioned according to IVAROLD
CP                         O: partitioned according to IVARNEW
CP   NORIG     I    INT    total number of variables (including fixed)
CP   IVAROLD   O    INT    information about old partitioning
CP                            i = 1..M      XORIG(IVAROLD(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAROLD(i)) independent
CP                            Note: fixed variables do not occur in IVAROLD
CP   IVARNEW   O    INT    information about old partitioning
CP                            i = 1..M      XORIG(IVARNEW(i)) dependent
CP                            i = (M+1)..N  XORIG(IVARNEW(i)) independent
CP                            Note: fixed variables do not occur in IVARNEW
CP   xtmp      W    DP
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
      integer NORIG
      double precision X(N)
      integer IVAROLD(N)
      integer IVARNEW(N)
      double precision xtmp(NORIG)
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
      do i = 1, N
         xtmp(IVAROLD(i)) = X(i)
      enddo
      do i = 1, N
         X(i) = xtmp(IVARNEW(i))
      enddo

      return
      end
