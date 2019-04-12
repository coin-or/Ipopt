C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine GET_SIGMA(NLB, NUB, S_L, S_U, V_L, V_U,
     1                     SIGMA_L, SIGMA_U)
C
C*******************************************************************************
C
C    $Id: get_sigma.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Compute primal and dual steps.
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
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   S_L       I    DP     slack variables for lower bounds
CP   S_U       I    DP     slack variables for upper bounds
CP   V_L       I    DP     dual variables for lower bounds
CP   V_U       I    DP     dual variables for upper bounds
CP   SIGMA_L   O    DP     primal-dual Hessian of lower bound barrier term
CP                            (NLB diagonal elements only)
CP   SIGMA_U   O    DP     primal-dual Hessian of upper bound barrier term
CP                            (NUB diagonal elements only)
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
      include 'IPOPT.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NLB
      integer NUB
      double precision S_L(NLB)
      double precision S_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision SIGMA_L(NLB)
      double precision SIGMA_U(NUB)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i

      double precision MPEC_ETA
      common /MPEC/ MPEC_ETA
      save /MPEC/
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( QMPEC_TRIGGER.ge.1 ) then
         do i = 1, NLB
            SIGMA_L(i) = V_L(i)/(S_L(i) + MPEC_ETA*V_L(i))
         enddo
         do i = 1, NUB
            SIGMA_U(i) = V_U(i)/(S_U(i) + MPEC_ETA*V_U(i))
         enddo
      else
         do i = 1, NLB
            SIGMA_L(i) = V_L(i)/S_L(i)
         enddo
         do i = 1, NUB
            SIGMA_U(i) = V_U(i)/S_U(i)
         enddo
      endif

      return
      end
