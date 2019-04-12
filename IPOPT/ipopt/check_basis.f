C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      logical function CHECK_BASIS(ITER, M, C, PY, CONDC)
C
C*******************************************************************************
C
C    $Id: check_basis.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Test for conditioning of basis matrix
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
CP   ITER      I    INT    iteration counter
CP   M         I    INT    number of constraints
CP   C         I    DP     values of constraints
CP   PY        I    DP     range space step
CP   CONDC     O    DP     estimated condition number of C
CP                            (only for output; is stored internally)
CP   CHECK_BASIS O  LOG    new basis required if true
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
CS    DNRM2
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
      integer ITER
      integer M
      double precision C(M)
      double precision PY(M)
      double precision CONDC
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      double precision RLARGE, RSMALL
      save             RLARGE, RSMALL

      double precision DNRM2

      double precision pynrm, ratio, cnrm, quot
      character*80 line(2)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      if( ITER .eq. 0 ) then
         RLARGE = 1.d-20
         RSMALL = 1.d20

         CHECK_BASIS = .false.
         CONDC = 0.d0
      endif

      pynrm = DNRM2(M, PY, 1)
      cnrm  = DNRM2(M, C , 1)

      if( pynrm .gt. 0 ) then
         ratio = cnrm/pynrm
      else
         ratio = RLARGE
      endif

      if( ITER.gt.0 .and. pynrm.gt.0 ) then
         RLARGE = dmax1(ratio, RLARGE)
         RSMALL = dmin1(ratio, RSMALL)

         quot = RLARGE/RSMALL
         write(line,700) ITER, quot
 700     format(/,' In Iter ',i5,
     1        ' Check Basis with RLARGE/RSMALL = ',
     1        d11.3)
         call C_OUT(1,2,2,line)

         CONDC = RLARGE/RSMALL
         CHECK_BASIS = (CONDC .gt. QMAXCOND)

         if( CHECK_BASIS ) then
            RLARGE = 1.d-20
            RSMALL = 1.d20
         endif

      else

         CHECK_BASIS = .false.

      endif

      return
      end
