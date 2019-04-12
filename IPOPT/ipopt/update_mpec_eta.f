C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine UPDATE_MPEC_ETA(N, M, X, C, LAM, NLB, V_L, NUB, 
     1                           V_U, MU, RESTO)
C*******************************************************************************
C
C    
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Update parameter used in modification of step calculation
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
CA    
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
CP   N         I    INT    number of variables (without fixed)
CP   M         I    INT    number of constraints
CP   X         I    DP     actual iterate (reordered without fixed vars:
CP                             first M entries belong to dependent
CP                             variables, remaining to independent variables)
CP   C         I    DP     values of constraints at X
CP   LAM       I    DP     current values of multipliers
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   ILB       I    INT    indices of lower bounds
CP   V_L       I    DP     dual variables for lower bounds
CP                            (e.g. BNDS_L(i) is bound for X(ILB(i)) )
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. BNDS_U(i) is bound for X(IUB(i)) )
CP   V_U       I    DP     dual variables for upper bounds
CP   MU        I    DP     barrier parameter
CP   RESTO     I    INT    if not 0, we are in restoration phase and set
CP                            MPEC_ETA to zero
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
CS     IDAMAX
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
      integer M
      double precision X(N)
      double precision C(M)
      double precision LAM(M)
      integer NLB
      double precision V_L(NLB)
      integer NUB
      double precision V_U(NUB)
      double precision MU
      integer RESTO
C
C*******************************************************************************
C
C                              Local Variables
C
C*******************************************************************************
C
      integer IDAMAX
      integer i
      double precision scale, LASTMU, CURRMU
      save scale, LASTMU, CURRMU

      character*128 line

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
 10   if( QMPEC_TRIGGER.eq.1 .and. RESTO.le.0 ) then
         if( MU.gt.QMPEC_THRESH ) then
            MPEC_ETA = 0.d0
            LASTMU = MU
            CURRMU = MU
         else
            scale = 0.d0
            if( M.gt.0 ) then
               i = IDAMAX(M,C,1)
               scale = dabs(C(i))
            endif
            if( NLB.gt.0 ) then
               i = IDAMAX(NLB,V_L,1)
               scale = dmax1(scale,V_L(i))
            endif
            if( NUB.gt.0 ) then
               i = IDAMAX(NUB,V_U,1)
               scale = dmax1(scale,V_U(i))
            endif
            scale = scale + 1
            if( CURRMU.ne.MU ) then
               LASTMU = CURRMU
               CURRMU = MU
            endif
            MPEC_ETA = QMPEC_ETAFACT*LASTMU/scale
         endif
         if( QCNR.gt.0 .and. QPRINT.ge.4 ) then
            write(line,8000) MPEC_ETA
 8000       format('MPEC_ETA set to ',d10.4)
            call C_OUT(1,0,1,line)
         endif
      elseif( QMPEC_TRIGGER.gt.1 ) then
         call C_OUT(2,0,1,'Nothing implemented for QMPEC_TRIGGER.gt.1')
         call C_OUT(2,0,1,'Proceeding with QMPEC_TRIGGER.eq.1')
         QMPEC_TRIGGER = 1
         goto 10
      else
         MPEC_ETA = 0.d0
      endif

      return
      end
