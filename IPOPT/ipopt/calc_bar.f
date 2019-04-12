C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      double precision function CALC_BAR(NLB, NUB, S_L, S_U,
     1                                   V_L, V_U, MU)
C
C*******************************************************************************
C
C     $Id: calc_bar.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Calculate barrier part of primal-dual merit function
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
CP   V_L       I    DP     dual variables for lower bounds (only for QMERIT=+-1)
CP   V_U       I    DP     dual variables for upper bounds (only for QMERIT=+-1)
CP   MU        I    DP     barrier parameter
CP   CALC_BAR  O    DP     barrier part of primal-dual merit function
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
CS    DDOT
CS    C_OUT
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
      double precision MU
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision DDOT, dlog

      integer i
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( abs(QMERIT).eq.1 ) then
         if( QALPHA.ne.0 ) then
            write(line,*) 'calc_bar: Need same step size for QMERIT = ',
     1           QMERIT
            call C_OUT(2,0,1,line)
            stop
         endif

         CALC_BAR = DDOT(NLB, S_L, 1, V_L, 1) +
     1              DDOT(NUB, S_U, 1, V_U, 1)

         do i = 1, NLB
            CALC_BAR = CALC_BAR - MU*dlog((S_L(i)**2)*V_L(i))
CCC            CALC_BAR = CALC_BAR - MU*(2*dlog(S_L(i))+dlog(V_L(i)))
         enddo

         do i = 1, NUB
            CALC_BAR = CALC_BAR - MU*dlog((S_U(i)**2)*V_U(i))
CCC            CALC_BAR = CALC_BAR - MU*(2*dlog(S_U(i))+dlog(V_U(i)))
         enddo
      elseif( abs(QMERIT).eq.2 .or. QMERIT.eq.3 .or. abs(QMERIT).eq.4
     1        .or. abs(QMERIT).eq.5 ) then
         CALC_BAR = 0.d0
         do i = 1, NLB
            CALC_BAR = CALC_BAR - MU*dlog(S_L(i))
         enddo

         do i = 1, NUB
            CALC_BAR = CALC_BAR - MU*dlog(S_U(i))
         enddo
      elseif( QMERIT.eq.0 ) then
         CALC_BAR = 0d0
      else
         write(line,*) 'Invalid choice of QMERIT = ',QMERIT
         call C_OUT(2,0,1,line)
         stop
      endif

      return
      end
