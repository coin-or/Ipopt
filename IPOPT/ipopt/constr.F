C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine CONSTR(TASK, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1     NORIG, XORIG, CSCALE, VIN, VOUT, IVEC1, IVEC2,
     2     LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: constr.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Dirty business of dealing with the constraints...
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
CP   TASK      I    INT    =0: Provide maximum for work space needs
CP                                LRSTORE, LISTORE, LRW, LIW
CP                                (first call only,
CP                                 i.e. flag for initialization)
CP                             Also, give NZORIG into IVEC(1) (needed for full
CP                                 version only)
CP                         =1: Determine partitioning (might also happen, when
CP                                basis becomes singular?)
CP                             Determine IVAR (ignoring fixed variables), so
CP                                that the columns in A with the numbers
CP                                IVAR(NIND+1..N) are linear independent
CP                                in other words:
CP                                i=1..M       XORIG(IVAR(i)) is dep.
CP                                i=X+1..N     XORIG(IVAR(i)) is indep.
CP                         =2: Compute values of constraints C at XORIG (store
CP                                in VOUT)
CP                         =3: Compute CC*VOUT = VIN, where CC are the columns
CP                                of A corresponding to the dependent variables.
CP                                This is for IVEC1(1) = 1. For IVEC1(1) = 0
CP                                solve CC^{T}*VOUT = VIN
CCP                         =3: Evaluate Jacobian at XORIG and solve for VOUT
CCP                                from CC*VOUT = VIN, where CC are the columns
CCP                                of A corresponding to the dependent variables.
CCP                                Also, compute CinvN = CC^{-1}NN and store it
CCP                                (NN columns of A corresponding to indep.
CCP                                 vars.)
CCP                                for later calls with TASK = 4,5,6 (can use RS
CCP                                and IS for this.)
CP                         =4: Compute VOUT = CinvN^T * VIN
CP                         =5: Compute VOUT = CinvN   * VIN
CP                         =6: Compute VOUT = CinvN' * diag(VIN) * CinvN
CP                             (Note, VOUT is expected in 'packed form',
CP                              i.e. only the upper diagonal elements
CP                              are to be stored:
CP                               VOUT(1)  =   row 1  col 1
CP                               VOUT(2)  =   row 1  col 2
CP                               VOUT(3)  =   row 2  col 2
CP                               VOUT(4)  =   row 1  col 3
CP                               VOUT(5)  =   row 2  col 3
CP                               VOUT(6)  =   row 3  col 3  etc.)
CP                         =7: Copy rows of CinvN into VOUT:
CP                             VOUT(i + NIND*(j-1)) = CinvN(IVEC1(i),j)
CP                             (only for IVEC1(i) > 0 !)
CP                                i = 1,...,NIND
CP                         =8: Compute VOUT = A*VIN (Here VOUT is ordered
CP                                as the partitioned X, not XORIG; i.e.
CP                                first all dependent, then all independent
CP                                variables)
CP                         =9: Compute VOUT = A'*VIN (Here VIN is ordered
CP                                as the partitioned X, not XORIG; i.e.
CP                                first all dependent, then all independent
CP                                variables)
CP                        =10: Evaluate Jacobian, ordered accoring to IVAR
CP                                (first all dependent vars)
CP                                VIN(1): number of NZA as reals
CP                                VIN(2): number of NZC as reals
CP                                           (I know it's bad)
CP                                VOUT: nonzero elements of A
CP                                IVEC1: ACON (row indices)
CP                                IVEC2: AVAR (column indices)
CP                        =11: Copy current values of Jacobian into internal
CP                             storage space (only possible if initialized
CP                             with QQUASI < 0 or QQUASI >5!) for later call
CP                             with TASK = 12.
CP                        =12: Compute VOUT = AOLD'*VIN (Here VIN is ordered
CP                                as the partitioned X, not XORIG; i.e.
CP                                first all dependent, then all independent
CP                                variables).  AOLD is the one from the latest
CP                                call with TASK = 11.  Note that this will
CP                                be garbage if IVAR changed in between!
CP                        =13: Evaluate Jacobian (don't store!) at XORIG and
CP                                compute VOUT = Anew * VIN.
CP                                Note: VOUT is ordered as XORIG and does
CP                                include entries corresponding to fixed
CP                                variables.  Don't scale accorind to CSCALE!
CP                        =14: Compute the maxnorm of each constraint gradient
CP                                Result is stored in VOUT(M).
CP                        =15: Return the number of nonzero elements of Jacobian
CP                                (returned in IVEC1(1))
CP                        =20: Compute product VOUT of the (weighted) Hessians
CP                                and VIN.  The weights are given by CSCALE.
CP                                (It is assumed here, that no further scaling
CP                                 is necessary...)
CP                                Here, VIN and VOUT are ordered like XORIG
CP   ITER      I    INT    iteration counter
CP                            (it is assumed that only one evaluation of A has
CP                             to be done per iteration, i.e. this can be
CP                             use to figure out, if A has to be reevaluated.
CP                             It is also usedto see if CC has already been
CP                             factorized in this iteration.)
CP   N         I    INT    number of not-fixed variables
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints = number of dependent vars
CP   IVAR     I/O   INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: usually, IVAR is mapping from partition to
CP                                  XORIG (as described in the last 2 lines).
CP                                  However, on return for TASK = 1, it is only
CP                                  the map from the partition to X(!!!) without
CP                                  taking care of fixed variables.  Those are
CP                                  taken care of in 'PARTITION'.
CP                            Output for TASK = 1, Input for other TASKS
CP   NFIX      I    INT    number of fixed variables
CP                            (only for TASK = 1,3)
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP                            (only for TASK = 1,3)
CP   NORIG     I    INT    total number of variables in problem statement
CP                            (incl. fixed variables)
CP   XORIG     I    DP     actual iterate
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   CSCALE    I    DP     scaling factors for constraints C.  All output has
CP                            to be scaled accoring to CSCALE.
CP                            CSCALE is a scalar, unless QSCALE = 2, then CSCALE
CP                            is a vector with scaling factors for the
CP                            individual constraints
CP                            Also, for TASK = 20 values of weights for Hessians
CP   VIN       I    DP     see TASK
CP   VOUT      O    DP     see TASK
CP   IVEC1    I/O   INT    TASK = 7: Info about what part of CinvN is needed
CP                               10: ACON for A
CP   IVEC2    I/O   INT    TASK =10: AVAR for A
CP   LRS      I/O   INT    length of RS (Output only for TASK = 0)
CP   RS       I/O   DP     can be used to store DP variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LIS      I/O   INT    length of IS (Output only for TASK = 0)
CP   IS       I/O   INT    can be used to store INT variables between calls;
CP                            this array is not touched from the calling
CP                            program
CP   LRW      I/O   INT    length of RW (Output only for TASK = 0)
CP                            (This program has to check, if LRW is indeed
CP                             large enough!)
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW      I/O   INT    length of IW (Output only for TASK = 0)
CP                            (This program has to check, if LIW is indeed
CP                             large enough!)
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective function
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
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
CS    MA28_CALL
CS    MA48_CALL
CS    EVAL_A
CS    EVAL_C
CS    CONSTR_NEW_A
CS    CONSTR_ORDER_A
CS    CONSTR_BACKS
CS    C_OUT
CS    DCOPY
CS    DSCAL
CS    DGEMV
CS    DSPR
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

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer TASK
      integer ITER
      integer NORIG
      integer N
      integer NIND
      integer M
      integer IVAR(*)
      integer NFIX
      integer IFIX(NFIX)
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      double precision VIN(*)
      double precision VOUT(*)
      integer IVEC1(*)
      integer IVEC2(*)
      integer LRS
      double precision RS(LRS)
      integer LIS
      integer IS(LIS)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer NZORIG, LRS_MA48, LIS_MA48, LRW_MA48, LIW_MA48
      integer P_A, P_RSMA48, P_CINVN, P_ACON, P_AVAR, P_ISMA48
      integer P_AOLD
      integer NZC, NZA
      save    NZORIG, LRS_MA48, LIS_MA48, LRW_MA48, LIW_MA48
      save    P_A, P_RSMA48, P_CINVN, P_ACON, P_AVAR, P_ISMA48
      save    P_AOLD
      save    NZC, NZA

C flags:
      integer LASTITER_A, LASTITER_C
      save    LASTITER_A, LASTITER_C

      logical MA48INIT
      save    MA48INIT

      integer QCNFLAG
C      parameter( QCNFLAG = 1 )

      integer p_iwend, p_rwend, p_is, p_rs, p_x, p_d
      integer p_isig, nsig, ia, nzn, p_itmp, p_n2a, p_ncol, p_scinvn
      integer p_ivarold, p_ivartmp, lfix, p_xcor, p_islacks, p_ixo2x
      integer p_atmp, p_avar2, p_acon2
      integer idummy, k, i, j, lsig, row, col, icor, nztrash, nzctrash
      double precision dummy, scal_val
      logical get_cinvn, ex_slacks

      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      QCNFLAG = 2
      p_iwend = 0
      p_rwend = 0
      IERR = 0

      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1     1200, 1300, 1400, 1500, 1600)
     1     TASK + 1

      goto (2000) TASK - 19

      write(line,*) 'Invalid flag for TASK in CONSTR = ',TASK
      call C_OUT(2,0,1,line)
      stop

 100  continue
C-------------------------------------------------------------------------------
C     Start: Initialize NZORIG, LRS_MA48, LIS_MA48, LRW_MA48, LIW_MA48
C            and pointers for A, AVAR, ACON, RSMA28, ISMA28
C-------------------------------------------------------------------------------

C     Get NZORIG, LRS_MA28, LIS_MA28, LRW_MA28, LIW_MA28
      call EV_A(0, NORIG, XORIG, NZORIG, dummy, idummy, idummy,
     1     DAT, IDAT)
      if( QFULL.eq.0 .or. QLAMBDA.eq.1 .or.
     1     (QQUASI.ne.0.and.abs(QQUASI).lt.6) .or.
     1     abs(QERROR).eq.1 .or. QLAMBDA.eq.0 ) then
C
C     TODO: In those cases no partition is needed and MA28 or MA48 should
C     never be called!
C
#ifndef HAVE_MA48
         call MA28_CALL(
#else
         call MA48_CALL(
#endif
     1        0, N, M, NZORIG, dummy, idummy, idummy, 0, 0,
     1        dummy, idummy, LRS_MA48, dummy, LIS_MA48, idummy,
     2        LRW_MA48, dummy, LIW_MA48, idummy, IERR)
         MA48INIT= .true.
      else
         LRS_MA48 = 0
         LIS_MA48 = 0
         LRW_MA48 = 0
         LIW_MA48 = 0
         MA48INIT= .false.
      endif
      P_A      = 0
      P_RSMA48 = P_A      + NZORIG
      if( QCNFLAG.eq.1 ) then
         P_CINVN  = P_RSMA48 + LRS_MA48
         LRS      = P_CINVN  + M*NIND
      else
         P_CINVN  = -100
         LRS      = P_RSMA48 + LRS_MA48
      endif
      if( QQUASI.lt.0 .or. QQUASI.ge.6 ) then
         P_AOLD = LRS
         LRS    = P_AOLD + NZORIG
      endif

      P_ACON   = 0
      P_AVAR   = P_ACON   + NZORIG
      P_ISMA48 = P_AVAR   + NZORIG
      LIS      = P_ISMA48 + LIS_MA48
C
C     Determine need of work space (not very tight yet...)
C

C     TASK = 1
      LRW = LRW_MA48+LRS_MA48
      LIW = max(NORIG,LIW_MA48+LIS_MA48+2*N+NORIG)
C     TASK = 3
      LRW = max(LRW, LRW_MA48+NZORIG)
      LIW = max(LIW, LIW_MA48+2*NZORIG)
C
      IVEC1(1) = NZORIG
C
C     Initialize memory of when last evaluation of Jacobian was done
C
      LASTITER_A = -1
      LASTITER_C = -1

C-------------------------------------------------------------------------------
C     End:   Initializing work space and computing LRSTORE, LISTORE, LRW, LIW
C-------------------------------------------------------------------------------
      goto 9999

 200  continue
C-------------------------------------------------------------------------------
C     Start: Determine IVAR (partitioning into dependent and independent
C            variables from constraints (PARTITION SOMEWHERE ELSE!) )
C-------------------------------------------------------------------------------

C
C     If there is no old evaluation of A, evaluate at new point, get rid of
C     fixed variable columns, and get partition.
C
CTODO:  Current implementation of the 'else' case is wrong, since
C       the result of the swapping sorting is different when started from
C       A which has been sorted before.  For now, just reevaluate Jacobian.
C       But it might be nicer to avoid this extra evaluation.
      if( ITER.ne.LASTITER_A .or. .true. ) then   ! need to do something
C      if( ITER.ne.LASTITER_A ) then
C
         p_islacks = p_iwend
         p_ixo2x   = p_islacks + N
         p_iwend   = p_ixo2x + NORIG
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         call CONSTR_NEW_A(1, NZORIG, NORIG, XORIG, M, NIND,
     1        IW(p_ixo2x+1), NFIX, IFIX,
     1        CSCALE, NZA, RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1),
     2        idummy, LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C,
     5        EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
            goto 9999
         endif
C
C     Determine which of the colums of A belong to slack variables and
C     scale them up, so that MA28 will choose them into basis
C
         inquire(file = 'SLACKS.DAT', exist = ex_slacks)
         if( ex_slacks ) then
            scal_val = 1d10
            do i = 1, N
               IW(p_islacks+i) = 0
            enddo
            open(10,file='SLACKS.DAT',status='unknown')
 205        read(10,'(i16)',end=210) j
            k = IW(p_ixo2x+j)
            if( k.ne.0 ) then
               IW(p_islacks+k) = 1
            endif
            goto 205
 210        continue
            close(10)
            p_rwend = p_ixo2x
            do i = 1, NZA
               j = IS(P_AVAR+i)
               k = IW(p_islacks+j)
               if( k.ne.0 ) then
                  RS(P_A+i) = RS(P_A+i)*scal_val
               endif
            enddo
         endif
C
C     Call MA28 or MA48 to do factorization and obtain partitioning in IVAR
C
         p_rs    = p_rwend
         p_rwend = p_rs + LRS_MA48
         p_is    = p_iwend
         p_iwend = p_is + LIS_MA48
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         elseif( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         if( .not.MA48INIT ) then
            call C_OUT(2,0,1,
     1           'constr: BUG (partition) MA48/28 not initialized!')
            IERR = 95
            goto 9999
         endif
#ifndef HAVE_MA48
         call MA28_CALL(
#else
         call MA48_CALL(
#endif
     1        1, N, M, NZA, RS(P_A+1), IS(P_ACON+1),
     1        IS(P_AVAR+1), 0, 0, dummy, IVAR,
     2        LRS_MA48, RW(p_rs+1), ! Don't need to store RS
     3        LIS_MA48, IW(p_is+1), ! Don't need to store IS
     4        LRW-p_rwend, RW(p_rwend+1),
     5        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.lt.0 ) then
            write(line,*)
#ifndef HAVE_MA48
     1          'Constr (partition): Warning in MA28_CALL: IERR = ',IERR
#else
     1          'Constr (partition): Warning in MA48_CALL: IERR = ',IERR
#endif
            call C_OUT(2,0,1,line)
            IERR = 0
         elseif( IERR.ne.0 ) then
            write(line,*)
#ifndef HAVE_MA48
     1           'Constr (Partition): Error in MA28_CALL: IERR = ', IERR
#else
     1           'Constr (Partition): Error in MA48_CALL: IERR = ', IERR
#endif
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         p_rwend = p_rs
         p_iwend = p_is
C
C     Undo scaling
C
         if( ex_slacks ) then
            do i = 1, NZA
               j = IS(P_AVAR+i)
               k = IW(p_islacks+j)
               if( k.ne.0 ) then
                  RS(P_A+i) = RS(P_A+i)/scal_val
               endif
            enddo
         endif
C
C     Finally, reorder A, so that C are the first entries
C
         call CONSTR_ORDER_A(NZA, N, M, NIND, IVAR,
     1        0, NZA, RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC,
     2        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*) 'CONSTR-1: CONSTR_ORDER_A returns IERR =',
     1           IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif

      else
C
C     If there is an old evaluation of A, use this to get partition (MA28),
C     correct wrong entries in AVAR, and correct IVAR
C
         write(*,*) 'This part  of CONSTR is not working!'
         STOP
C
C     First, store old IVAR, but take dependent vars first
C
         p_ivarold = p_iwend
         p_iwend   = p_ivarold + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         do i = 1, M
            IW(p_ivarold+i) = IVAR(i+NIND)
         enddo
         do i = 1, NIND
            IW(p_ivarold+i+M) = IVAR(i)
         enddo
C
C     Now take out fixed vars
C
         if( NFIX.gt.0 ) then
            p_xcor = p_iwend  ! Here store by what entries have to be corrected
            p_iwend = p_xcor + NORIG
            if( p_iwend.gt.LIW ) then
               IERR = 99
               goto 9999
            endif
            lfix = 1
            icor = 0
            do i = 1, NORIG
               if( i.eq.IFIX(lfix) ) then
                  icor = icor + 1
                  lfix = lfix + 1
                  if( lfix.gt.NFIX ) then
                     do j = i+1, NORIG
                        IW(p_xcor+j) = icor
                     enddo
                     goto 220
                  endif
               else
                  IW(p_xcor+i) = icor
               endif
            enddo
 220        continue
            do i = 1, N
               j = IW(p_ivarold+i)
               IW(p_ivarold+i) = j - IW(p_xcor+j)
            enddo
            p_iwend = p_xcor
         endif
C
C     Store new partition in ivartmp (since it referres to wrong ACON)
C
         p_ivartmp = p_iwend
         p_iwend   = p_ivartmp + N
c$$$C
c$$$C     Now get partition from MA28
c$$$C
c$$$         p_rs    = p_rwend
c$$$         p_rwend = p_rs + LRS_MA28
c$$$         p_is    = p_iwend
c$$$         p_iwend = p_is + LIS_MA28
c$$$         if( p_rwend.gt.LRW ) then
c$$$            IERR = 98
c$$$            goto 9999
c$$$         elseif( p_iwend.gt.LIW ) then
c$$$            IERR = 99
c$$$            goto 9999
c$$$         endif
c$$$         call MA28_CALL(1, N, M, NZA, RS(P_A+1), IS(P_ACON+1),
c$$$     1        IS(P_AVAR+1), 0, 0, dummy, IW(p_ivartmp+1),
c$$$     2        LRS_MA28, RW(p_rs+1), ! Don't need to store RS
c$$$     3        LIS_MA28, IW(p_is+1), ! Don't need to store IS
c$$$     4        LRW-p_rwend, RW(p_rwend+1),
c$$$     5        LIW-p_iwend, IW(p_iwend+1), IERR)
c$$$         if( IERR.lt.0 ) then
c$$$            write(line,*)
c$$$     1          'Constr (partition): Warning in MA28_CALL: IERR = ',IERR
c$$$            call C_OUT(2,0,1,line)
c$$$            IERR = 0
c$$$         elseif( IERR.ne.0 ) then
c$$$            write(line,*)
c$$$     1           'Constr (Partition): Error in MA28_CALL: IERR = ', IERR
c$$$            call C_OUT(2,0,1,line)
c$$$            IERR = 10
c$$$            goto 9999
c$$$         endif
c$$$         p_rwend = p_rs
c$$$         p_iwend = p_is
C
C     Now, reorder A, so that ACON is corrected and the first entries are C
C
         call CONSTR_ORDER_A(NZA, N, M, NIND, IW(p_ivartmp+1),
     1        0, NZA, RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC,
     2        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*) 'CONSTR-1: CONSTR_ORDER_A returns IERR =',
     1           IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
C
C     Finally, correct IVAR
C
         do i = 1, N
            IVAR(i) = IW(p_ivarold+IW(p_ivartmp+i))
         enddo

      endif
C
C     Reset Flags
C
      LASTITER_A = ITER
      LASTITER_C = -1

C-------------------------------------------------------------------------------
C     End:   Determine IVAR (partitioning into dependent and independent
C            variables
C-------------------------------------------------------------------------------
      goto 9999

 300  continue
C-------------------------------------------------------------------------------
C     Start: Evaluate constraints C
C-------------------------------------------------------------------------------

      call EV_C(NORIG, XORIG, M, VOUT, DAT, IDAT)
      if( QSCALE.ge.2 ) then
         do i = 1, M
            VOUT(i) = VOUT(i)*CSCALE(i)
         enddo
      elseif( CSCALE(1).ne.1d0 ) then
         call DSCAL(M, CSCALE(1), VOUT, 1)
      endif

C-------------------------------------------------------------------------------
C     End:   Evaluate constraints C
C-------------------------------------------------------------------------------
      goto 9999

 400  continue
C-------------------------------------------------------------------------------
C     Start: Evaluate Jacobian A and solve for specified vectors
C-------------------------------------------------------------------------------

C
C     First determine whether CinvN will have to be recomputed as well:
C
      if( QCNFLAG.eq.1 .and. (LASTITER_C.ne.ITER) ) then
         get_cinvn = .true.
      else
         get_cinvn = .false.
      endif
C
C     Solve system with C
C
      call DCOPY(M, VIN, 1, VOUT, 1)
      call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1     NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2     RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, IVEC1(1), VOUT,
     2     LRS_MA48, LIS_MA48, RS(P_RSMA48+1), IS(P_ISMA48+1),
     3     LRW_MA48, LIW_MA48, LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend,
     4     IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.gt.0 ) then
         call C_OUT(2,0,1,'CONSTR-3: Error from CONSTR_BACKS')
         goto 9999
      endif
C
C     If QCNFLAG = 1 and CinvN has not been computed before in this iteration
C     compute CinvN
C
      if( get_cinvn ) then
C
C     Store N in RS(P_CINVN+1)
C
         call DCOPY( NIND*M, 0.d0, 0, RS(P_CINVN+1), 1)
         do i = NZC+1, NZA
            col = IS(P_AVAR+i) - M
            row = IS(P_ACON+i)
            RS(P_CINVN + M*(col-1) + row) = RS(P_A+i)
         enddo
C
C     Get CinvN from MA48:
C
         if( .not.MA48INIT ) then
            call C_OUT(2,0,1,
     1           'constr: BUG (CinvN) MA48/28 not initilized!')
            IERR =95
            goto 9999
         endif
#ifndef HAVE_MA48
         call MA28_CALL(
#else
         call MA48_CALL(
#endif
     1        3, N, M, NZC, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), 1, NIND, RS(P_CINVN+1), idummy,
     3        LRS_MA48, RS(P_RSMA48+1),
     4        LIS_MA48, IS(P_ISMA48+1),
     5        LRW-p_rwend, RW(p_rwend+1),
     6        LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.lt.0 ) then
            write(line,*)
#ifndef HAVE_MA48
     1           'Constr (CinvN): Warning in MA28_CALL: IERR = ', IERR
#else
     1           'Constr (CinvN): Warning in MA48_CALL: IERR = ', IERR
#endif
            call C_OUT(2,0,1,line)
            IERR = 0
         elseif( IERR.ne.0 ) then
            write(line,*)
#ifndef HAVE_MA48
     1           'Constr (CinvN): Error in MA28_CALL: IERR = ', IERR
#else
     1           'Constr (CinvN): Error in MA48_CALL: IERR = ', IERR
#endif
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif

C-------------------------------------------------------------------------------
C     End:   Evaluate Jacobian A and solve for specified vectors
C-------------------------------------------------------------------------------
      goto 9999

 500  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT (with CinvN transpose)  (using last evaluation of A)
C-------------------------------------------------------------------------------

      if( QCNFLAG.eq.1 ) then
         if( LASTITER_C.ne.ITER ) then
            write(line,510) LASTITER_C, ITER
 510        format('Constr (4): Error with CinvN:',/,
     1             '            LASTITER_C = ',i5,', but ITER = ',i5)
            call C_OUT(2,0,2,line)
            stop
         endif
         call DGEMV('T', M, NIND, 1.d0, RS(P_CINVN+1), M, VIN, 1,
     1        0.d0, VOUT, 1)
      else
         p_x     = p_rwend
         p_rwend = p_x + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     Solve C^T x = vin
C
         call DCOPY(M, VIN, 1, RW(p_x+1), 1)
         call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1        NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2        RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 0, RW(p_x+1),
     2        LRS_MA48, LIS_MA48, RS(P_RSMA48+1), IS(P_ISMA48+1),
     3        LRW_MA48, LIW_MA48, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR-4: Error from CONSTR_BACKS')
            goto 9999
         endif
C
C     Compute vout = N^T x
C
         call DCOPY(NIND, 0d0, 0, VOUT, 1)
         do i = NZC+1, NZA
            row = IS(P_AVAR+i) - M
            col = IS(P_ACON+i)
            VOUT(row) = VOUT(row) + RS(P_A+i)*RW(p_x+col)
         enddo
         p_rwend = p_x
      endif

C-------------------------------------------------------------------------------
C     End:   Compute VOUT (with CinvN transpose)  (using last evaluation of A)
C-------------------------------------------------------------------------------

      goto 9999

 600  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT (with CinvN)  (using last evaluation of A)
C-------------------------------------------------------------------------------

      if( QCNFLAG.eq.1 ) then
         if( LASTITER_C.ne.ITER ) then
            write(line,610) LASTITER_C, ITER
 610        format('Constr (6): Error with CinvN:',/,
     1             '            LASTITER_C = ',i5,', but ITER = ',i5)
            call C_OUT(2,0,2,line)
            IERR = 4
            goto 9999
         endif
         call DGEMV('N', M, NIND, 1.d0, RS(P_CINVN+1), M, VIN, 1,
     1        0.d0, VOUT, 1)
      else

C
C     Compute vout = N vin
C
         call DCOPY(M, 0d0, 0, VOUT, 1)
         do i = NZC+1, NZA
            col = IS(P_AVAR+i) - M
            row = IS(P_ACON+i)
            VOUT(row) = VOUT(row) + RS(P_A+i)*VIN(col)
         enddo
C
C     Solve C vout = vout
C
         call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1        NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2        RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 1, VOUT,
     2        LRS_MA48, LIS_MA48, RS(P_RSMA48+1), IS(P_ISMA48+1),
     3        LRW_MA48, LIW_MA48, LRW-p_rwend, RW(p_rwend+1),
     4        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR-4: Error from CONSTR_BACKS')
            goto 9999
         endif
      endif

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Compute VOUT (with CinvN)  (using last evaluation of A)
C-------------------------------------------------------------------------------

 700  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT = CinvN' * diag(VIN) * CinvN
C            (in packed form!)
C-------------------------------------------------------------------------------

C
C     Determine indices in VIN that are non-zero
C      (ISIG contains indices of nonzeros):
C
      p_isig   = p_iwend
      p_iwend = p_isig   + M
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      nsig = 0
      do i = 1, M
         if( VIN(i).ne.0.d0 ) then
            nsig = nsig + 1
            IW(p_isig+nsig) = i
         endif
      enddo
      if( nsig.eq.0 ) then     ! Nothing to do!
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
         goto 9999
      endif
      p_iwend = p_isig + nsig

      if( QCNFLAG.eq.1 ) then
C
C     Compute BB (in packed form!)
C
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
C
C     Do nsig rank-one updates:
C
         do i = 1, nsig
            lsig = IW(p_isig+i)
            call DSPR('U', NIND, VIN(lsig), RS(P_CINVN+lsig), M, VOUT)
         enddo

      elseif( QCNFLAG.eq.0 .or. (QCNFLAG.eq.2.and.nsig.le.NIND) ) then

C
C     Use only adjoint factorizations and rank-one updates
C

C
C     Initialize matrix (could also just add these to B, need to change get_bb)
C
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
C
C     Do rank one-updates for each nonzero in VIN (Sigma_d)
C
         p_x     = p_rwend
         p_d     = p_x + M
         p_rwend = p_d + NIND
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         do i = 1, nsig
C
C     Set x = e_i
C
            lsig = IW(p_isig+i)
            call DCOPY(M, 0.d0, 0, RW(p_x+1), 1)
            RW(p_x+lsig) = 1.d0
C
C     Solve C^T x = x
C
            call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1           NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2           RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 0,
     2           RW(p_x+1), LRS_MA48, LIS_MA48, RS(P_RSMA48+1),
     3           IS(P_ISMA48+1), LRW_MA48, LIW_MA48, LRW-p_rwend,
     4           RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2           EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     5           DAT, IDAT)
            if( IERR.gt.0 ) then
               call C_OUT(2,0,1,'CONSTR-6: Error from CONSTR_BACKS')
               goto 9999
            endif
C
C     Compute d = N^T x
C
            call DCOPY(NIND, 0d0, 0, RW(p_d+1), 1)
            do j = NZC+1, NZA
               row = IS(P_AVAR+j) - M
               col = IS(P_ACON+j)
               RW(p_d+row) = RW(p_d+row) + RS(P_A+j)*RW(p_x+col)
            enddo
C
C     Do the rank-one update
C
            call DSPR('U', NIND, VIN(lsig), RW(p_d+1), 1, VOUT)

         enddo

         p_rwend = p_x

      elseif( QCNFLAG.eq.2 ) then
C
C     Initialize matrix (could also just add these to B, need to change get_bb)
C
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
C
C     In this case, we compute CinvN by sending all columns of N to MA48,
C     but we will store only those belonging to ISIGS in SCINVN ...
C
         p_scinvn = p_rwend
         p_x      = p_scinvn + nsig*NIND
         p_rwend  = p_x + M
         p_ncol   = p_iwend
         p_n2a    = p_ncol + NIND  ! n2a gives relation between sorted N and A
         p_iwend  = p_n2a + NZA-NZC
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         elseif( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
C
C     In order to access columns of N fast, order N into n, ncol, nvar
C
         do i = 1, NIND
            IW(p_ncol+i) = 0
         enddo
C     count nonzeros in each column
         do i = NZC+1, NZA
            col = IS(P_AVAR+i)-M
            IW(p_ncol+col) = IW(p_ncol+col) + 1
         enddo
C     now sort
         p_itmp  = p_iwend
         p_iwend = p_itmp + NIND
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         IW(p_itmp+1) = 0
         do i = 2, NIND
            IW(p_itmp+i) = IW(p_ncol+i-1) + IW(p_itmp+i-1)
         enddo
         do ia = NZC+1, NZA
            col = IS(P_AVAR+ia)-M
            IW(p_itmp+col) = IW(p_itmp+col) + 1
            nzn = IW(p_itmp+col)
            IW(p_n2a+nzn) = ia
         enddo
         p_iwend = p_itmp
C
         nzn = 0
         do i = 1, NIND
C
C     Store i=th column of N in x
C
            call DCOPY(M, 0d0, 0, RW(p_x+1), 1)
            do j = nzn+1, nzn+IW(p_ncol+i)
               ia = IW(p_n2a+j)
               RW(p_x+IS(P_ACON+ia)) = RS(P_A+ia)
            enddo
            nzn = nzn + IW(p_ncol+i)
C
C     Call MA48 to solve C x = x
C
            call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1           NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2           RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 1,
     1           RW(p_x+1),
     2           LRS_MA48, LIS_MA48, RS(P_RSMA48+1), IS(P_ISMA48+1),
     3           LRW_MA48, LIW_MA48, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               call C_OUT(2,0,1,'CONSTR-4: Error from CONSTR_BACKS')
               goto 9999
            endif
C
C     Copy out significant entries of SCINVN
C
            do j = 1, nsig
               lsig = IW(p_isig+j)
               RW(p_scinvn+j+(i-1)*nsig) = RW(p_x+lsig)
            enddo
         enddo
C
C     Do nsig rank-one updates:
C
         do i = 1, nsig
            lsig = IW(p_isig+i)
            call DSPR('U', NIND, VIN(lsig), RW(p_scinvn+i), nsig, VOUT)
         enddo
         p_rwend = p_scinvn
         p_iwend = p_ncol

      else
         call C_OUT(2,0,1,'CONSTR-6: QCNFLAG invalid!')
      endif

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute VOUT = CinvN' * diag(VIN) * CinvN
C            (in packed form!)
C-------------------------------------------------------------------------------

 800  continue
C-------------------------------------------------------------------------------
C     Start: Copy part of CinvN into VOUT
C-------------------------------------------------------------------------------

      if( QCNFLAG.eq.1 ) then
         if( LASTITER_C.ne.ITER ) then
            call C_OUT(2,0,1,'CONSTR-7: LASTITER_C wrong!')
            stop
         endif
C
C     Copy stuff directly out of CINVN
C
         do i = 1, NIND
            k = IVEC1(i)
            if( k.gt.0 ) then
               call DCOPY(NIND, RS(P_CINVN+k), M, VOUT(i), NIND)
            endif
         enddo
      elseif( QCNFLAG.eq.0 .or. QCNFLAG.eq.2 ) then
C
C     Compute requested rows of CinvN by adjoint factorization
C
         p_x     = p_rwend
         p_rwend = p_x + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         do i = 1, NIND
            k = IVEC1(i)
            if( k.gt.0 ) then
C
C     Set x = e_k
C
               call DCOPY(M, 0.d0, 0, RW(p_x+1), 1)
               RW(p_x+k) = 1.d0
C
C     Solve C^T x = x
C
               call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1              NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE,
     2              NZA, RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 0,
     2              RW(p_x+1), LRS_MA48, LIS_MA48, RS(P_RSMA48+1),
     3              IS(P_ISMA48+1), LRW_MA48, LIW_MA48, LRW-p_rwend,
     4              RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR,
     2              EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     5              EV_HCV, DAT, IDAT)
               if( IERR.gt.0 ) then
                  call C_OUT(2,0,1,'CONSTR-7: Error from CONSTR_BACKS')
                  goto 9999
               endif
C
C     Compute d = N^T x
C
               call DCOPY(NIND, 0d0, 0, VOUT(i), NIND)
               do j = NZC+1, NZA
                  row = IS(P_AVAR+j) - M
                  col = IS(P_ACON+j)
                  VOUT(i+(row-1)*NIND) = VOUT(i+(row-1)*NIND) +
     1                 RS(P_A+j)*RW(p_x+col)
               enddo
            endif
         enddo

      elseif( QCNFLAG.eq.2 ) then
C
C     Compute all columns of CinvN after each other and store what is needed in
C     VOUT...
C

C
C     In this case, we compute CinvN by sending all columns of N to MA48,
C     but we will store only those needed in VOUT
C
         p_x      = p_rwend
         p_rwend  = p_x + M
         p_ncol   = p_iwend
         p_n2a    = p_ncol + NIND  ! n2a gives relation between sorted N and A
         p_iwend  = p_n2a + NZA-NZC
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         elseif( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
C
C     In order to access columns of N fast, order N into n, ncol, nvar
C
         do i = 1, NIND
            IW(p_ncol+i) = 0
         enddo
C     count nonzeros in each column
         do i = NZC+1, NZA
            col = IS(P_AVAR+i)-M
            IW(p_ncol+col) = IW(p_ncol+col) + 1
         enddo
C     now sort
         p_itmp  = p_iwend
         p_iwend = p_itmp + NIND
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         IW(p_itmp+1) = 0
         do i = 2, NIND
            IW(p_itmp+i) = IW(p_ncol+i-1) + IW(p_itmp+i-1)
         enddo
         do ia = NZC+1, NZA
            col = IS(P_AVAR+ia)-M
            IW(p_itmp+col) = IW(p_itmp+col) + 1
            nzn = IW(p_itmp+col)
            IW(p_n2a+nzn) = ia
         enddo
         p_iwend = p_itmp
C
         nzn = 0
         do i = 1, NIND
C
C     Store i=th column of N in x
C
            call DCOPY(M, 0d0, 0, RW(p_x+1), 1)
            do j = nzn+1, nzn+IW(p_ncol+i)
               ia = IW(p_n2a+j)
               RW(p_x+IS(P_ACON+ia)) = RS(P_A+ia)
            enddo
            nzn = nzn + IW(p_ncol+i)
C
C     Call MA48 to solve C x = x
C
            call CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1           NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA,
     2           RS(P_A+1), IS(P_ACON+1), IS(P_AVAR+1), NZC, 1,
     1           RW(p_x+1),
     2           LRS_MA48, LIS_MA48, RS(P_RSMA48+1), IS(P_ISMA48+1),
     3           LRW_MA48, LIW_MA48, LRW-p_rwend, RW(p_rwend+1),
     4           LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
            if( IERR.gt.0 ) then
               call C_OUT(2,0,1,'CONSTR-4: Error from CONSTR_BACKS')
               goto 9999
            endif
C
C     Copy out significant entries into VOUT
C
            do j = 1, NIND
               k = IVEC1(j)
               if( k.gt.0 ) then
                  VOUT(j+(i-1)*NIND) = RW(p_x+k)
               endif
            enddo
         enddo

         p_rwend = p_x
         p_iwend = p_ncol

      else

         call C_OUT(2,0,1,'CONSTR-7: Invalid QCNFLAG!')
         IERR = 4
         goto 9999
      endif
      goto 9999

C-------------------------------------------------------------------------------
C     End:   Copy part of CinvN into VOUT
C-------------------------------------------------------------------------------

 900  continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A*Vin
C-------------------------------------------------------------------------------

      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR-8: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
      call DCOPY(N, 0.d0, 0, VOUT, 1)
C
C     A * Vin
C
      do i = 1, NZA
         row = IS(P_AVAR+i)
         col = IS(P_ACON+i)
         VOUT(row) = VOUT(row) + RS(P_A+i)*VIN(col)
      enddo

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute Vout = A*Vin
C-------------------------------------------------------------------------------

 1000 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A'*Vin
C-------------------------------------------------------------------------------

      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR-9: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
      call DCOPY(N, 0.d0, 0, VOUT, 1)
C
C     A^T * Vin
C
      do i = 1, NZA
         col = IS(P_AVAR+i)
         row = IS(P_ACON+i)
         VOUT(row) = VOUT(row) + RS(P_A+i)*VIN(col)
      enddo

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute Vout = A'*Vin
C-------------------------------------------------------------------------------

 1100 continue
C-------------------------------------------------------------------------------
C     Start: Compute Jacobian
C-------------------------------------------------------------------------------

C
C     If necessary, compute A (and store in RS, IS)
C
      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     Copy into output variables
C
      call DCOPY( NZA, RS(P_A+1), 1, VOUT, 1)
      do i = 1, NZA
         IVEC1(i) = IS(P_ACON+i)
      enddo
      do i = 1, NZA
         IVEC2(i) = IS(P_AVAR+i)
      enddo
      vin(1) = dble(NZA)
      vin(2) = dble(NZC)

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Compute Jacobian
C-------------------------------------------------------------------------------

 1200 continue
C-------------------------------------------------------------------------------
C     Start: Copy current Jacobian into storage space
C-------------------------------------------------------------------------------

      if( QQUASI.ge.0 .and. QQUASI.lt.6 ) then
         IERR = 4
         goto 9999
      endif
C
C     If necessary, compute A (and store in RS, IS)
C
      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     Copy into storage space
C
      call DCOPY( NZA, RS(P_A+1), 1, RS(P_AOLD+1), 1)

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Copy current Jacobian into store space
C-------------------------------------------------------------------------------

 1300 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = AOLD*Vin
C-------------------------------------------------------------------------------

      if( QQUASI.ge.0 .and. QQUASI.lt.6 ) then
         IERR = 4
         goto 9999
      endif
      call DCOPY(N, 0.d0, 0, VOUT, 1)
C
C     AOLD * Vin
C
      do i = 1, NZA
         row = IS(P_AVAR+i)
         col = IS(P_ACON+i)
         VOUT(row) = VOUT(row) + RS(P_AOLD+i)*VIN(col)
      enddo

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute Vout = AOLD*Vin
C-------------------------------------------------------------------------------

 1400 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A(at XORIG) * Vin
C-------------------------------------------------------------------------------

      p_atmp  = p_rwend
      p_rwend = p_atmp + NZORIG
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      p_avar2 = p_iwend
      p_acon2 = p_avar2 + NZORIG
      p_iwend = p_acon2 + NZORIG
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
C
C     Evaluate Jacobian at XORIG
C
      call CONSTR_NEW_A(0, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1     NFIX, IFIX, CSCALE, nztrash, RW(p_atmp+1), IW(p_acon2+1),
     1     IW(p_avar2+1), nzctrash, LIW-p_iwend, IW(p_iwend+1), IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1     DAT, IDAT)
      if( IERR.gt.0 ) then
         call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
         goto 9999
      endif
C
C     Do the product
C
      call DCOPY(NORIG, 0.d0, 0, VOUT, 1)
      do i = 1, NZORIG
         row = IW(p_avar2+i)
         col = IW(p_acon2+i)
         VOUT(row) = VOUT(row) + RW(p_atmp+i)*VIN(col)
      enddo

      p_rwend = p_atmp
      p_iwend = p_avar2

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute Vout = A(at XORIG) * Vin
C-------------------------------------------------------------------------------

 1500 continue
C-------------------------------------------------------------------------------
C     Start: Compute maxnorm of each constraint gradient
C-------------------------------------------------------------------------------

C
C     If necessary, compute A (and store in RS, IS)
C
      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     Copy into storage space
C
      call DCOPY(M, 0.d0, 0, VOUT, 1)
      do i = 1, NZORIG
         row = IS(P_ACON+i)
         VOUT(row) = max(VOUT(row), abs(RS(P_A+i)))
      enddo

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute maxnorm of each constraint gradient
C-------------------------------------------------------------------------------

1600  continue
C-------------------------------------------------------------------------------
C     Start: Return number of nonzero elements in Jacobian
C-------------------------------------------------------------------------------

C
C     If necessary, compute A (and store in RS, IS)
C
      if( LASTITER_A.ne.ITER ) then
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, IFIX, CSCALE, NZA, RS(P_A+1), IS(P_ACON+1),
     2        IS(P_AVAR+1), NZC, LIW-p_iwend, IW(p_iwend+1), IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3        DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR: Error from CONSTR_NEW_A!')
            goto 9999
         endif
         LASTITER_A = ITER
      endif
      IVEC1(1) = NZA

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Return number of nonzero elements in Jacobian
C-------------------------------------------------------------------------------

 2000 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = (weighted Hessian of constraints) * Vin
C-------------------------------------------------------------------------------

      call EV_HCV(0, NORIG, XORIG, M, CSCALE, VIN, VOUT, DAT, IDAT)

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Compute Vout = (weighted Hessian of constraints) * Vin
C-------------------------------------------------------------------------------


 9999 continue
      return
      end


C*******************************************************************************
C
      subroutine CONSTR_NEW_A(TASK, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1     NFIX, IFIX, CSCALE, NZ, A, ACON, AVAR, NZC, LIW, IW, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     2     DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: constr.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Evaluate Jacobian, eliminate colums for fixed variables, and
CT    re-order (first C, then N)
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
CP   TASK      I    INT    =0: Evaluate A and reorder and eliminate nothing
CP                             (don't even scale rows!)
CP                         =1: Evaluate A, eliminate columns corresponding
CP                             to fixed entries.  Don't re-order A otherwise.
CP                         =2: Evaluate A, eliminate columns corresponding
CP                             to fixed entries and order A, so that first
CP                             NZC entries correspond to basis matrix C, rest
CP                             to N.
CP   NZORIG    I    INT    nubmer of non-zeros in original A (including fixed)
CP   NORIG     I    INT    number of all variables (including fixed ones)
CP   XORIG     I    DP     point, where A has to be evaluated
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   M         I    INT    number of constraints = number of dependent vars
CP                            (only for TASK=2)
CP   NIND      I    INT    number of independent variables
CP                            (only for TASK=2)
CP   IVAR     I/O   INT    TASK = 1: Output
CP                            XORIG(i) corresponds to X(IVAR(i)) if X not
CP                            reordered (only fixed vars taken out)
CP                         TASK = 2: Input
CP                            information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP                            (only for TASK=1)
CP   CSCALE    I    DP     scaling factors for cosntraints C.  All output has
CP                            to be scaled accoring to CSCALE.
CP                            CSCALE is a scalar, unless QSCALE >= 2, then
CP                            CSCALE is a vector with scaling factors for the
CP                            individual constraints.  If = 3, then 
CP                            CSCALE(M+1...M+N) constains the scaling values
CP                            for the variables
CP   NZ        O    INT    number of non-zeros in A after eliminating columns
CP   A         O    INT    non-zero entries of A after removing columns for
CP                            for fixed variables.
CP                            TASK = 1:
CP                               numbering of columns as in original Jacobian
CP                               but without fixed columns
CP                            TASK = 2:
CP                               numbering accoring to partition;
CP                               first M columns correspond to dependent
CP                               variables, rest to independent variables
CP   ACON      O    INT    row indices of A
CP   AVAR      O    INT    column indices of A
CP   NZC       O    INT    number of nonzeros in basis matrix C
CP                            (only for TASK=2)
CP   LIW       I    INT    length of IW
CP   IW        W    INT    integer work space
CP   IERR      O    INT    IERR not 0: ERROR
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective function
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
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
CS    DSCAL
CS    EVAL_A
CS    CONSTR_ORDER_A
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
      integer TASK
      integer NZORIG
      integer NORIG
      double precision XORIG(NORIG)
      integer M
      integer NIND
      integer NFIX
      integer IVAR(*)
      integer IFIX(NFIX)
      double precision CSCALE(*)
      integer NZ
      double precision A(NZORIG)
      integer ACON(NZORIG)
      integer AVAR(NZORIG)
      integer NZC
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, k, p_iwend, lfix, lx
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_iwend = 0
      IERR = 0
C
C     Evaluate Jacobian
C
      call EV_A(1, NORIG, XORIG, NZORIG, A, ACON, AVAR, DAT, IDAT)
C
C     For TASK = 1 (unpartitioned)
C
C     Get rid of colums that belong to fixed variables
C     and determine number of non-zero element in this matrix
C     (assume that IFIX is ordered in increasing order!)
C
      if( TASK.eq.1 ) then
         if( NFIX.gt.0 ) then
C     Find out, which of the XORIG are not fixed
            lfix = 1
            lx = 1
            do i = 1, NORIG
               if( lfix.le.NFIX ) then
                  if( i.eq.IFIX(lfix) ) then
                     IVAR(i) = 0
                     lfix = lfix + 1
                     goto 100
                  endif
               endif
               IVAR(i) = lx
               lx = lx + 1
 100           continue
            enddo
            NZ = 0
            do i = 1, NZORIG
               k = IVAR(AVAR(i))
               if( k.gt.0 ) then
                  NZ = NZ + 1
                  AVAR(NZ) = k
                  ACON(NZ) = ACON(i)
                  A   (NZ) = A(i)
               endif
            enddo
         else
            do i = 1, NORIG
               IVAR(i) = i
            enddo
            NZ = NZORIG
         endif
C
C     For TASK = 2 (IVAR given)
C
C     Get rid of colums that belong to fixed variables and re-order
C     A so that the first NZC (output) entries belong to C and the rest to N.
C
      elseif( TASK.eq.2 ) then

         call CONSTR_ORDER_A(NZORIG, NORIG, M, NIND, IVAR, NFIX, NZ,
     1        A, ACON, AVAR, NZC, LIW-p_iwend, IW(p_iwend+1), IERR)
         if( IERR.ne.0 ) then
            write(line,*)
     1           'CONSTR_NEW_A: CONSTR_ORDER_A returns IERR =',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
C
      elseif( TASK.ne.0 ) then
C
C     Invalid TASK
         call C_OUT(2,0,1,'Wrong TASK in constr_new_a')
         STOP
      endif
C
C     scale rows accoring to CSCALE
C
      if( TASK.ne.0 ) then
         if( QSCALE.eq.2 ) then
            do i = 1, NZ
               A(i) = A(i) * CSCALE(ACON(i))
            enddo
         elseif( QSCALE.ge.3 ) then
            do i = 1, NZ
               A(i) = A(i) * CSCALE(ACON(i)) * CSCALE(M+AVAR(i))
            enddo
         elseif( CSCALE(1).ne.1d0 ) then
            call DSCAL(NZ, CSCALE(1), A, 1)
         endif
      endif
cCDELETEME
c      do i = 1, NZ
c         write(QCNR,*) 'A(',ACON(i),',',AVAR(i),') = ',A(i)
c      enddo

 9999 continue
      end

C*******************************************************************************
C
      subroutine CONSTR_BACKS(ITER, LASTITER_A, LASTITER_C,
     1     NZORIG, NORIG, XORIG, M, NIND, NFIX, IVAR, CSCALE, NZA, A,
     2     ACON, AVAR, NZC, TRANSP, V, LRS_MA48, LIS_MA48,
     3     RS_MA48, IS_MA48, LRW_MA48, LIW_MA48, LRW, RW, LIW, IW,
     4     IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2     EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: constr.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Do a backsolve involving C or C^T
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
CP   ITER      I    INT    actual iterate
CP   LASTITER_A I/O INT    I: last iterate, when A has been evaluated
CP                         O: ITER
CP   LASTITER_C I/O INT    I: last iterate, when C has been factorized
CP                         O: ITER
CP   NZORIG    I    INT    number of non-zeros in original A (including fixed)
CP   NORIG     I    INT    number of all variables (including fixed ones)
CP   XORIG     I    DP     point, where A has to be evaluated
CP                            XORIG is ordered in ORIGINAL order (i.e. not
CP                            partitioned into independent and dependent
CP                            variables)
CP   M         I    INT    number of constraints = number of dependent vars
CP                            (only for TASK=2)
CP   NIND      I    INT    number of independent variables
CP                            (only for TASK=2)
CP   NFIX      I    INT    number of fixed variables
CP   IVAR     I/O   INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP                            Output for TASK = 1, Input for TASK = 3
CP   CSCALE    I    DP     scaling factors for cosntraints C.  All output has
CP                            to be scaled accoring to CSCALE.
CP                            CSCALE is a scalar, unless QSCALE = 2, then CSCALE
CP                            is a vector with scaling factors for the
CP                            individual constraints
CP   NZA       O    INT    number of non-zeros in A after eliminating columns
CP   A         O    INT    non-zero entries of A after removing columns for
CP                            for fixed variables.
CP                               numbering accoring to partition;
CP                               first M columns correspond to dependent
CP                               variables, rest to independent variables
CP   ACON      O    INT    row indices of A
CP   AVAR      O    INT    column indices of A
CP   NZC       O    INT    number of nonzeros in basis matrix C
CP                            (only for TASK=2)
CP   TRANSP    I    INT    =1: solve with C, otherwise with C^T
CP   V        I/O   DP     I: rhs  O: solution
CP   LRS_MA48  I    INT
CP   LIS_MA48  I    INT
CP   RS_MA48   S    DP
CP   IS_MA48   S    INT
CP   LRW_MA48  I    INT
CP   LIW_MA48  I    INT
CP   LRW       I    INT    length of RW
CP   RW        W    INT    real work space
CP   LIW       I    INT    length of IW
CP   IW        W    INT    integer work space
CP   IERR      O    INT    IERR not 0: ERROR
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective function
CP   EV_A      I    EXT    Subroutine for Jacobian
CP   EV_H      I    EXT    Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT    Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT    Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT    Subroutine for constraint Hessian-vector products
CP   DAT       P    DP     privat DP data for evaluation routines
CP   IDAT      P    INT    privat INT data for evaluation routines
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
CS    CONSTR_NEW_A
CS    MA48_CALL
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
      integer ITER
      integer LASTITER_A
      integer LASTITER_C
      integer NZORIG
      integer NORIG
      double precision XORIG(NORIG)
      integer M
      integer NIND
      integer NFIX
      integer IVAR(NORIG-NFIX)
      double precision CSCALE(*)
      integer NZA
      double precision A(NZORIG)
      integer ACON(NZORIG)
      integer AVAR(NZORIG)
      integer NZC
      integer TRANSP
      double precision V(M)
      integer LRS_MA48
      integer LIS_MA48
      double precision RS_MA48(LRS_MA48)
      integer IS_MA48(LIS_MA48)
      integer LRW_MA48
      integer LIW_MA48
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer IERR
      external EV_F
      external EV_C
      external EV_G
      external EV_A
      external EV_H
      external EV_HLV
      external EV_HOV
      external EV_HCV
      double precision DAT(*)
      integer IDAT(*)
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer p_rwend, p_iwend, job, idummy
      character*80 line
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_rwend = 0
      p_iwend = 0
      IERR = 0
C
C     Need to evaluate Jacobian only when not moved
C
      if( ITER.ne.LASTITER_A ) then
         LASTITER_A = ITER
         call CONSTR_NEW_A(2, NZORIG, NORIG, XORIG, M, NIND, IVAR,
     1        NFIX, idummy, CSCALE, NZA, A, ACON, AVAR, NZC, LIW, IW,
     2        IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2        EV_HCV, DAT, IDAT)
         if( IERR.gt.0 ) then
            call C_OUT(2,0,1,'CONSTR_BACKS: Error from CONSTR_NEW_A!')
            goto 9999
         endif
      endif
C
C     Get VOUT from MA48
C
      IF( LASTITER_C.ne.ITER ) then
         job = 2
      else
         job = 3
      endif
#ifndef HAVE_MA48
      call MA28_CALL(
#else
      call MA48_CALL(
#endif
     1               job, M+NIND, M, NZC, A, ACON,
     1               AVAR, TRANSP, 1, V, idummy,
     2               LRS_MA48, RS_MA48,
     3               LIS_MA48, IS_MA48,
     4               LRW-p_rwend, RW(p_rwend+1),
     5               LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.lt.0 ) then
         write(line,*)
     1      'CONSTR_BACKS: Warning in MA48_CALL: IERR = ', IERR
         call C_OUT(2,0,1,line)
         IERR = 0
      elseif( IERR.ne.0 ) then
         write(line,*)
     1      'CONSTR_BACKS: Error in MA48_CALL: IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
      LASTITER_C = ITER
C
 9999 continue
      end

C*******************************************************************************
C
      subroutine CONSTR_ORDER_A(NZORIG, NORIG, M, NIND, IVAR,
     1     NFIX, NZ, A, ACON, AVAR, NZC, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: constr.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Order A (first C, then N)
C
C-------------------------------------------------------------------------------
C                          Programm description
C-------------------------------------------------------------------------------
C
CB    Order A, so that elements for basis matrix C are first.
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
CP   NZORIG    I    INT    number of non-zeros in original A (including fixed)
CP   NORIG     I    INT    number of all variables (including fixed ones)
CP   M         I    INT    number of constraints = number of dependent vars
CP                            (only for TASK=2)
CP   NIND      I    INT    number of independent variables
CP                            (only for TASK=2)
CP   IVAR     I/O   INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NFIX      I    INT    number of fixed variables
CP   NZ        O    INT    number of non-zeros in A after eliminating columns
CP   A        I/O   INT    non-zero entries of A after removing columns for
CP                            for fixed variables.
CP                               numbering accoring to partition;
CP                               first M columns correspond to dependent
CP                               variables, rest to independent variables
CP                            Input: unordered  Output: ordered
CP   ACON     I/O   INT    row indices of A
CP                            Input: unordered  Output: ordered
CP   AVAR     I/O   INT    column indices of A
CP                            Input: unordered  Output: ordered
CP   NZC       O    INT    number of nonzeros in basis matrix C
CP                            (only for TASK=2)
CP   LIW       I    INT    length of IW
CP   IW        W    INT    integer work space
CP   IERR      O    INT    IERR not 0: ERROR
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
      integer NZORIG
      integer NORIG
      integer M
      integer NIND
      integer NFIX
      integer IVAR(NORIG-NFIX)
      integer NZ
      double precision A(NZORIG)
      integer ACON(NZORIG)
      integer AVAR(NZORIG)
      integer NZC
      integer LIW
      integer IW(LIW)
      integer IERR
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer i, k, p_ivar1, p_iwend
      integer lc1, lc2, ln1, ln2, acon_n, avar_n, ndiff
      double precision a_n
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
CTODO Don't to the work if not necessary!

      p_iwend = 0
      IERR = 0
C
C     Get rid of colums that belong to fixed variables and re-order
C     A so that the first NZC (output) entries belong to C and the rest to N.
C

C
C     First, compute the inverse of IVAR
C
      p_ivar1 = p_iwend
      p_iwend = p_ivar1 + NORIG
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      if( NFIX.gt.0 ) then
         do i = 1, NORIG
            IW(p_ivar1+i) = 0
         enddo
      endif
      do i = 1, NORIG-NFIX
         IW(p_ivar1+IVAR(i)) = i
      enddo
C
C     Now correct column indices of A, get rid of columns corresponding to
C     fixed variables, and reorder A
C
      lc1 = 1                   ! first free entry for C
      lc2 = 0                   ! last examined entry in A (upwards)
      ln1 = NZORIG              ! first free entry for N (counting backwards)
      ln2 = NZORIG+1            ! last examined entry in A (downwards)
C
C     Find next entry in from below that belongs to N
C
 110  continue
      lc2 = lc2 + 1
      if( lc2.eq.ln2 ) then
         goto 200
      endif
      k = IW(p_ivar1+AVAR(lc2))
      if( k.ne.0 ) then         ! not fixed variable
         if( k.le.M ) then   ! dependent variable
            A   (lc1) = A   (lc2)
            ACON(lc1) = ACON(lc2)
            AVAR(lc1) = k
            lc1 = lc1 + 1
            goto 110
         else                   ! This entry belongs to N
            a_n    = A(lc2)
            acon_n = ACON(lc2)
            avar_n = k
         endif
      else
         goto 110
      endif
C
C     Find next entry from above that belongs to C
C
 120  continue
      ln2 = ln2 - 1
      if( lc2.eq.ln2 ) then     ! visited all elements of A
         A   (ln1) = a_n
         ACON(ln1) = acon_n
         AVAR(ln1) = avar_n
         ln1 = ln1 - 1
         goto 200
      endif
      k = IW(p_ivar1+AVAR(ln2))
      if( k.ne.0 ) then         ! not fixed
         if( k.gt.M ) then   ! belongs to N
            A   (ln1) = A   (ln2)
            ACON(ln1) = ACON(ln2)
            AVAR(ln1) = k
            ln1 = ln1 - 1
            goto 120
         else                   ! This entry belongs to C: swap
            A   (lc1) = A   (ln2)
            ACON(lc1) = ACON(ln2)
            AVAR(lc1) = k
            lc1 = lc1 + 1
            A   (ln1) = a_n
            ACON(ln1) = acon_n
            AVAR(ln1) = avar_n
            ln1 = ln1 - 1
            goto 110
         endif
      else
         goto 120
      endif
C
C     Now close the gap caused by elimination of fixed columns
C
 200  continue
      ndiff = ln1 - lc1 + 1
      if( ndiff.gt.0 ) then
         do i = ln1+1, NZORIG
            A   (i-ndiff) = A   (i)
            ACON(i-ndiff) = ACON(i)
            AVAR(i-ndiff) = AVAR(i)
         enddo
      endif
C
C     Now determine number of nonzeros in C and in overall A
C
      NZC = lc1 - 1
      NZ = NZORIG - ndiff

 9999 continue
      end
