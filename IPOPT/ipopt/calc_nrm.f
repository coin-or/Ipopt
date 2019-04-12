C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      double precision FUNCTION CALC_NRM(N, VEC)
C
C*******************************************************************************
C
C    $Id: calc_nrm.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Calculate norm of vector (usually of constraints) - what particular
CT    norm depends on QCNRM parameter
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
CP   CALC_NRM  O    DP     norm of VEC -
CP                           QCNRM = 1: ell_1 norm
CP                           QCNRM = 2: ell_2 norm
CP                           QCNRM = 3: max norm
C
CP   N         I    INT    length of vector VEC
CP   VEC       I    DP     elements of vector whose norm is to calculated
CP   FILTER_PHI I   DP     filter entries - phi part
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
      integer N
      double precision VEC(N)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, IDAMAX
      double precision DNRM2, DASUM
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( QCNRM.eq.1 ) then
         CALC_NRM = DASUM(N, VEC, 1)
      elseif( QCNRM.eq.2 ) then
         CALC_NRM = DNRM2(N, VEC, 1)
      elseif( QCNRM.eq.3 ) then
         if( N.eq.0 ) then
            CALC_NRM = 0.d0
         else
            i = IDAMAX(N, VEC, 1)
            CALC_NRM = dabs(VEC(i))
         endif
      else
         call C_OUT(2,0,1,'calc_nrm: Invalid value of QCNRM. Abort.')
         stop
      endif

 9999 continue
      return
      end
