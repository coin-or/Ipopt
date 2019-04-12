C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C*******************************************************************************
C
      subroutine CONSTR(TASK, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1                  NORIG, XORIG, CSCALE, VIN, VOUT, IVEC1, IVEC2,
     2                  LRS, RS, LIS, IS, LRW, RW, LIW, IW, IERR)
C
C*******************************************************************************
C
C    $Id: constr.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Interface from IPOPT to DAE2NLP
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
CA    Andreas Waechter      10/10/00
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
CP                         =1: Determine IVAR, so
CP                                that the columns in A with the numbers
CP                                IVAR(NIND+1..N) are linear independent
CP                                in other words:
CP                                i=1..M       XORIG(IVAR(i)) is dep.
CP                                i=X+1..N     XORIG(IVAR(i)) is indep.
CP                             NOTE: This has to be called in any case before
CP                                doing anything with TASK > 1!
CP                                Also, no fixed variables are allowed
CP                         =2: Compute values of constraints C at XORIG (store
CP                                in VOUT)
CP                         =3: Compute CC*VOUT = VIN, where CC are the columns
CP                                of A corresponding to the dependent variables.
CP                                This is for IVEC1(1) = 1. For IVEC1(1) = 0
CP                                solve CC^{T}*VOUT = VIN
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
CP                             NOT YET IMLEMENTED!
CP                         =7: Copy rows of CinvN into VOUT:
CP                             VOUT(i + NIND*(j-1)) = CinvN(IVEC1(i),j)
CP                             (only for IVEC1(i) > 0 !)
CP                                i = 1,...,NIND
CP                             NOT YET IMLEMENTED!
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
CP                             with QQUASI < 0!) for later call with TASK = 12.
CP                        =12: Compute VOUT = AOLD'*VIN (Here VIN is ordered
CP                                as the partitioned X, not XORIG; i.e.
CP                                first all dependent, then all independent
CP                                variables).  AOLD is the one from the latest
CP                                call with TASK = 12.  Note that this will
CP                                be garbage if IVAR changed in between!
CP                        =13: Evaluate Jacobian (don't store!) at XORIG and
CP                                compute VOUT = Anew * VIN.
CP                                Note: VOUT is ordered as XORIG and does
CP                                include entries corresponding to fixed
CP                                variables.  Don't scale accorind to CSCALE!
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
CP   CSCALE    I    DP     scaling factors for cosntraints C.  All output has
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
CS    C_OUT
CS    DAXPY
CS    DCOPY
CS    DAE2NLP
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
      include 'DYNAUX.INC'      ! Here NZ, NCOL etc are defined
C !DEC$ ATTRIBUTES DLLEXPORT :: /PARAMS/, /DYNAUX/
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
C
C-------------------------------------------------------------------------------
C                            Local varibales
C-------------------------------------------------------------------------------
C
      integer IDAMAX
      double precision DNRM2
      integer NNZA, LRSD2N, LISD2N, LRWD2N, LIWD2N
      save    NNZA, LRSD2N, LISD2N, LRWD2N, LIWD2N
      integer P_RSD2N, P_ISD2N, P_A, P_AOLD
      save    P_RSD2N, P_ISD2N, P_A, P_AOLD
      integer LASTITER_A, LASTITER_C
      save    LASTITER_A, LASTITER_C

      integer p_iwend, p_rwend, p_rw, p_iw, p_tmp, p_xtmp, p_atmp
      integer p_isig, p_ei, p_tmp2, lsig
      integer iflag, idummy, i, j, nsig, lw, nnzac, nu_opt, np_opt
      double precision dummy, fillin, h, vnrm
      logical ex
      character*80 line

      logical QHESSFINDIFF
      parameter( QHESSFINDIFF = .false. )

C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      p_iwend = 0
      p_rwend = 0
      IERR = 0

      NU_OPT = NU - NU_PROF
      NP_OPT = NP - NP_FIX

      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1     1200, 1300, 1400)
     1     TASK + 1

      goto (2000) TASK - 19

      write(line,*) 'Invalid flag for TASK in CONSTR = ',TASK
      call C_OUT(2,0,1,line)
      stop

 100  continue
C-------------------------------------------------------------------------------
C     Start: Initialize NNZA, LASTITER's, LRS, LIS, LRW, LIW
C-------------------------------------------------------------------------------

C
C     Check parameters
C
      if( QSELBAS.ne.0 ) then
         IERR = 4
         goto 9999
      endif
C
C     We need to make sure that no variable is fixed
C
      if( N.ne.NORIG ) then
         call C_OUT(2,0,1,'CONSTR: N and NORIG different!')
         IERR = 633
         goto 9999
      endif
C
C     Get sizes of storage and work space
C
      inquire(file = 'FILLIN.DAT', exist = ex)
      if( ex ) then
         open(7, file='FILLIN.DAT', status='old', ERR = 9008)
         read(7,'(d24.16)', ERR = 9008, END = 9008) fillin
      else
         fillin  = 20d0         ! so much fill-in is expected
      endif
      iflag   = 1               ! fast factorzation (iflag=0) is to be used?

      call ADDCON_NNZ(NZ, NY, NU_OPT, nnzac)
      LRWD2N  = (NZ+NY)*(2*NZ+NY+NU_OPT+NP_OPT) + nnzac +
     1     NCOL*(NZ+NY+NU_OPT)
      p_rw    = p_rwend
      p_rwend = p_rw + LRWD2N
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      LIWD2N  = NZ + NCOL*(NZ+NY+NU_OPT) + NP_OPT + 2*nnzac
      p_iw    = p_iwend
      p_iwend = p_iw + LIWD2N
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif

      call DAE2NLP(0, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT,
     1     XORIG, iflag, idummy, NNZA, dummy, fillin, dummy, dummy,
     2     LRSD2N, dummy, LISD2N, idummy, LRWD2N, RW(p_rw+1),
     3     LIWD2N, IW(p_iw+1), IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-0 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

      P_A     = 0
      P_RSD2N = P_A     + NNZA
      LRS     = P_RSD2N + LRSD2N
      if( QQUASI.lt.0 ) then
         P_AOLD = LRS
         LRS    = P_AOLD + NNZA
      endif

      P_ISD2N = 0
      LIS     = P_ISD2N + LISD2N

      LRW = LRWD2N
      LIW = LIWD2N
C
      IVEC1(1) = 0              ! Hmmm... where is that needed?
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
C     No scaling is allowed here
C
      if( QSCALE.eq.2 .or. CSCALE(1).ne.1d0 ) then
         call C_OUT(2,0,1,'CONSTR: Can''t do scaling of constraints.')
         IERR = 4
         goto 9999
      endif
C
C     Check, if there is a file 'PARTITION.DAT'.  If so, read this to obtain
C     partition.  Otherwise use MA28 within DAE2NLP to get the partitioning.
C
      inquire(file = 'PARTITION.DAT', exist = ex)
      if( ex ) then
         iflag = 0
         open(7, file='PARTITION.DAT', status='old')
         read(7,'(i10)') IVAR(1)
         if( IVAR(1).gt.0 ) then
            do i = 2, NCOL*(NZ+NY+NU_OPT)
               read(7,'(i10)') IVAR(i)
            enddo
         else
            do i = 1, NCOL*(NZ+NY+NU_OPT)
               IVAR(i) = i
            enddo
         endif
      else
         if( CRIT_ELE.lt.1 .or. CRIT_ELE.gt.NELE ) then
            IERR = 634
            call C_OUT(2,0,1,'CONSTR: Invalid value of CRIT_ELE.')
            goto 9999
         endif
         iflag = CRIT_ELE
      endif

      call DAE2NLP(1, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG, iflag, IVAR, 0, dummy, dummy, dummy, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1), LRW-p_rwend,
     3     RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-1 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Reset Flags
C
      LASTITER_A = -1
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

      call DAE2NLP(2, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG, iflag, IVAR, 0, dummy, dummy, VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1), LRW-p_rwend,
     3     RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-2 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
CWEG
c      do i = 1, M
c         write(*,*) i,VOUT(I)
c      enddo
c      pause

C-------------------------------------------------------------------------------
C     End:   Evaluate constraints C
C-------------------------------------------------------------------------------
      goto 9999

 400  continue
C-------------------------------------------------------------------------------
C     Start: Solve CC * VOUT = VIN   or   CC^T * VOUT = VIN
C-------------------------------------------------------------------------------

C
C     If necessary, reevaluate Jacobian
C
      if( ITER.ne.LASTITER_A ) then
         call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT,
     1        XORIG, iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
CWEG
c         do i = 1, NNZA
c            write(*,*) i,RS(P_A+i)
c         enddo
c         pause
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-3 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     If necessary, factorize Jacobian
C
      if( ITER.ne.LASTITER_C ) then
         if( CRIT_ELE.lt.1 .or. CRIT_ELE.gt.NELE ) then
            IERR = 634
            call C_OUT(2,0,1,'CONSTR: Invalid value of CRIT_ELE.')
            goto 9999
         endif
         iflag = CRIT_ELE
         call DAE2NLP(8, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT,
     1        XORIG, iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-8 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         LASTITER_C = ITER
      endif
C
C     Solve system
C
      if( IVEC1(1).eq.0 ) then
         iflag = 10
      else
         iflag = 9
      endif
      call DAE2NLP(iflag, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT,
     1     XORIG, idummy, IVAR, NNZA, RS(P_A+1), VIN, VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-9/10 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

C-------------------------------------------------------------------------------
C     End:   Solve CC * VOUT = VIN   or   CC^T * VOUT = VIN
C-------------------------------------------------------------------------------
      goto 9999

 500  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT = Cinv^T * VIN
C-------------------------------------------------------------------------------

      if( LASTITER_C.ne.ITER .or. LASTITER_A.ne.ITER ) then
         write(line,510) LASTITER_A, LASTITER_C, ITER
 510     format('Constr (4): Error with CinvN:',/,
     1        'LASTITER_A = ',i5,' LASTITER_C = ',i5,', but ITER = ',i5)
         call C_OUT(2,0,2,line)
         stop
      endif
C
C     Compute tmp = C^{-T} * VIN
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DAE2NLP(10, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT,
     1     XORIG, idummy, IVAR, NNZA, RS(P_A+1), VIN, RW(p_tmp+1),
     2     dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-10 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Compute VOUT = N^T * tmp
C
CWEG      write(*,*) 'before DAE2NLP-7'
      call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     idummy, IVAR, NNZA, RS(P_A+1), RW(p_tmp+1), VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
CWEG      write(*,*) 'after DAE2NLP-7'
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-10 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      p_rwend = p_tmp

C-------------------------------------------------------------------------------
C     End:   Compute VOUT = Cinv^T * VIN
C-------------------------------------------------------------------------------
      goto 9999

 600  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT = Cinv * VIN
C-------------------------------------------------------------------------------

      if( LASTITER_C.ne.ITER .or. LASTITER_A.ne.ITER ) then
         write(line,610) LASTITER_A, LASTITER_C, ITER
 610     format('Constr (5): Error with CinvN:',/,
     1        'LASTITER_A = ',i5,' LASTITER_C = ',i5,', but ITER = ',i5)
         call C_OUT(2,0,2,line)
         stop
      endif
C
C     Compute tmp = N * vin
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DAE2NLP(6, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     idummy, IVAR, NNZA, RS(P_A+1), VIN, RW(p_tmp+1), dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-6 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Compute VOUT = C^{-1} * tmp
C
      call DAE2NLP(9, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     idummy, IVAR, NNZA, RS(P_A+1), RW(p_tmp+1), VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-9 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

C-------------------------------------------------------------------------------
C     End:   Compute VOUT = Cinv * VIN
C-------------------------------------------------------------------------------
      goto 9999

 700  continue
C-------------------------------------------------------------------------------
C     Start: Compute VOUT = CinvN' * diag(VIN) * CinvN
C            (in packed form!)
C-------------------------------------------------------------------------------

      if( LASTITER_C.ne.ITER .or. LASTITER_A.ne.ITER ) then
         write(line,710) LASTITER_A, LASTITER_C, ITER
 710     format('Constr (6): Error with CinvN:',/,
     1        'LASTITER_A = ',i5,' LASTITER_C = ',i5,', but ITER = ',i5)
         call C_OUT(2,0,2,line)
         stop
      endif
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
      inquire(file='expensive',exist=ex) ! this is for time comparison
      do i = 1, M
         if( VIN(i).ne.0.d0 .or. ex) then
            nsig = nsig + 1
            IW(p_isig+nsig) = i
         endif
      enddo
      if( nsig.eq.0 ) then     ! Nothing to do!
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
         goto 9999
      endif

      if( nsig.lt.2*NIND ) then
C
C     Use only adjoint factorizations and rank-one updates
C

C
C     Initialize matrix
C
         call DCOPY( (NIND*(NIND+1)/2), 0.d0, 0, VOUT, 1)
C
C     Do rank one-updates for each nonzero in VIN (Sigma_d)
C
         p_tmp   = p_rwend
         p_tmp2  = p_tmp + max(M,NIND) !p_x
         p_rwend = p_tmp2 + M !p_d
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         do i = 1, nsig
C
C     Set tmp = e_i
C
            lsig = IW(p_isig+i)
            call DCOPY(M, 0.d0, 0, RW(p_tmp+1), 1)
            RW(p_tmp+lsig) = 1.d0
C
C     Compute tmp2 = C^{-T} * tmp
C
            call DAE2NLP(10, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE,
     1           TI, ZINIT, XORIG, idummy, IVAR, NNZA, RS(P_A+1),
     1           RW(p_tmp+1), RW(p_tmp2+1),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-10 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Compute tmp = N^T * tmp2
C
            call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1           ZINIT, XORIG, idummy, IVAR, NNZA, RS(P_A+1),
     1           RW(p_tmp2+1), RW(p_tmp+1),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-7 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Do the rank-one update
C
            call DSPR('U', NIND, VIN(lsig), RW(p_tmp+1), 1, VOUT)

         enddo

         p_rwend = p_tmp

      else
         p_iwend = p_isig + nsig
C
C     Compute products with unit vectors
C
         lw      = 1
         p_ei    = p_rwend
         p_tmp   = p_ei    + NIND
         p_tmp2  = p_tmp   + M
         p_rwend = p_tmp2  + M
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
         call DCOPY(NIND, 0d0, 0, RW(p_ei+1), 1)
         do i = 1, NIND
            RW(p_ei+i) = 1d0
            if( i.gt.1 ) RW(p_ei+i-1) = 0d0
C
C     Compute tmp = N * ei
C
            call DAE2NLP(6, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1           ZINIT, XORIG, idummy, IVAR, NNZA, RS(P_A+1),
     1           RW(p_ei+1), RW(p_tmp+1),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-6 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Compute tmp2 = C^{-1} * tmp
C
            call DAE2NLP(9, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1           ZINIT, XORIG, idummy, IVAR, NNZA, RS(P_A+1),
     1           RW(p_tmp+1), RW(p_tmp2+1),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-9 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Compute tmp2 = VIN * tmp2
C
            do j = 1, M
               RW(p_tmp2+j) = RW(p_tmp2+j) * VIN(j)
            enddo
C
C     Compute tmp = C^{-T} * tmp2
C
            call DAE2NLP(10, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE,
     1           TI, ZINIT, XORIG, idummy, IVAR, NNZA, RS(P_A+1),
     1           RW(p_tmp2+1), RW(p_tmp+1),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-10 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
C
C     Compute W(column) = N^T * tmp
C
            call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1           ZINIT, XORIG,
     1           idummy, IVAR, NNZA, RS(P_A+1), RW(p_tmp+1), VOUT(lw),
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-10 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            lw = lw + i
         enddo
      endif
      goto 9999

C-------------------------------------------------------------------------------
C     End:   Compute VOUT = CinvN' * diag(VIN) * CinvN
C            (in packed form!)
C-------------------------------------------------------------------------------
      goto 9999

 800  continue
C-------------------------------------------------------------------------------
C     Start: Copy part of CinvN into VOUT
C-------------------------------------------------------------------------------

      call C_OUT(2,0,1,'CONSTR 7: Can''t get part of CinvN yet')
      IERR = 2346
      goto 9999


C-------------------------------------------------------------------------------
C     End:   Copy part of CinvN into VOUT
C-------------------------------------------------------------------------------
      goto 9999

 900  continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A*Vin
C-------------------------------------------------------------------------------

C
C     If necessary, reevaluate Jacobian
C
      if( ITER.ne.LASTITER_A ) then
         call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT, XORIG,
     1        iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-3 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     Compute VOUT(dep) = C^T * VIN
C
      call DAE2NLP(5, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_A+1), VIN, VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-5 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Compute VOUT(indep) = N^T * VIN
C
      call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_A+1), VIN, VOUT(M+1), dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-7 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

C-------------------------------------------------------------------------------
C     End:   Compute Vout = A*Vin
C-------------------------------------------------------------------------------
      goto 9999

 1000 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A'*Vin
C-------------------------------------------------------------------------------

C
C     If necessary, reevaluate Jacobian
C
      if( ITER.ne.LASTITER_A ) then
         call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT,
     1        XORIG, iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-3 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     tmp = C * VIN(dep)
C
      p_tmp   = p_rwend
      p_rwend = p_tmp + M
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      call DAE2NLP(4, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_A+1), VIN, RW(p_tmp+1), dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-4 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     VOUT = N * VIN
C
      call DAE2NLP(6, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_A+1), VIN(M+1), VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-6 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     VOUT = VOUT + tmp
C
      call DAXPY(M, 1d0, RW(p_tmp+1), 1, VOUT, 1)

C-------------------------------------------------------------------------------
C     End:   Compute Vout = A'*Vin
C-------------------------------------------------------------------------------
      goto 9999

 1100 continue
C-------------------------------------------------------------------------------
C     Start: Compute Jacobian
C-------------------------------------------------------------------------------

      call C_OUT(2,0,1,'CONSTR: Can''t export values of A.')
      IERR = 4
      goto 9999

C-------------------------------------------------------------------------------
C     End:   Compute Jacobian
C-------------------------------------------------------------------------------

 1200 continue
C-------------------------------------------------------------------------------
C     Start: Copy current Jacobian into storage space
C-------------------------------------------------------------------------------

      if( QQUASI.ge.0 ) then
         IERR = 4
         goto 9999
      endif
C
C     If necessary, reevaluate Jacobian
C
      if( ITER.ne.LASTITER_A ) then
         call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT,
     1        XORIG, iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-3 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         LASTITER_A = ITER
      endif
C
C     Copy into storage space
C
      call DCOPY(NNZA, RS(P_A+1), 1, RS(P_AOLD+1), 1)

      goto 9999

C-------------------------------------------------------------------------------
C     End:   Copy current Jacobian into store space
C-------------------------------------------------------------------------------

 1300 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = AOLD*Vin
C-------------------------------------------------------------------------------

      if( QQUASI.ge.0 ) then
         IERR = 4
         goto 9999
      endif
C
C     Compute VOUT(dep) = C^T * VIN
C
      call DAE2NLP(5, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_AOLD+1), VIN, VOUT, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-5 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Compute VOUT(indep) = N^T * VIN
C
      call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RS(P_AOLD+1), VIN, VOUT(M+1), dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-7 return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif

C-------------------------------------------------------------------------------
C     End:   Compute Vout = AOLD*Vin
C-------------------------------------------------------------------------------
      goto 9999

 1400 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = A(at XORIG) * Vin
C-------------------------------------------------------------------------------

      p_atmp  = p_rwend
      p_xtmp  = p_atmp + NNZA
      p_rwend = p_xtmp + NORIG
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
C
C     Evaluate Jacobian at XORIG
C
      call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RW(p_atmp+1), dummy, dummy, dummy,
     2     LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-3tmp return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     Do the product
C
      call DAE2NLP(5, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RW(p_atmp+1), VIN, RW(p_xtmp+1),
     2     dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-5tmp return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1     ZINIT, XORIG,
     1     iflag, IVAR, NNZA, RW(p_atmp+1), VIN, RW(p_xtmp+M+1),
     2     dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4     IERR)
      if( IERR.gt.0 ) then
         write(line,*) 'CONSTR: DAE2NLP-7tmp return IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      do i = 1, NORIG
         VOUT(IVAR(i)) = RW(p_xtmp+i)
      enddo
      p_rwend = p_atmp

      goto 9999
C-------------------------------------------------------------------------------
C     End:   Compute Vout = A(at XORIG) * Vin
C-------------------------------------------------------------------------------

 2000 continue
C-------------------------------------------------------------------------------
C     Start: Compute Vout = (weighted Hessian of constraints) * Vin
C-------------------------------------------------------------------------------

CWEG
c      do i = 1, M
c         write(*,*) i,CSCALE(i)
c      enddo
c      stop
      if( QHESSFINDIFF ) then
CWEG
	WRITE(*,*) 'THIS IS FINDIFF IN CONSTR!!!!'
C
C     Compute approximation of Hessian - Vector product by finite difference
C     of Jacobians
C
C
C     If necessary, reevaluate Jacobian
C
         if( ITER.ne.LASTITER_A ) then
            call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1           ZINIT,
     1           XORIG, iflag, IVAR, NNZA, RS(P_A+1), dummy, dummy,
     2           dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3           LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4           IERR)
            if( IERR.gt.0 ) then
               write(line,*) 'CONSTR: DAE2NLP-3 return IERR = ',IERR
               call C_OUT(2,0,1,line)
               goto 9999
            endif
            LASTITER_A = ITER
         endif
         p_atmp  = p_rwend
         p_xtmp  = p_atmp + NNZA
         p_rwend = p_xtmp + NORIG
         if( p_rwend.gt.LRW ) then
            IERR = 98
            goto 9999
         endif
C
C     Determine step size
C
CTODO         h = 1d-8
         i = IDAMAX(NORIG, VIN, 1)
         vnrm = dabs(VIN(i))
c         vnrm = DNRM2(NORIG, VIN, 1)
         if( vnrm.eq.0.d0 ) then
            call DCOPY(NORIG, 0.d0, 0, VOUT, 1)
            goto 9999
         else
            h = 1d-8/vnrm
         endif
         h = 1d-6
C
C     Compute Jacobian at X + h * VIN
C
         call DCOPY(NORIG, XORIG, 1, RW(p_xtmp+1), 1)
         call DAXPY(NORIG, h, VIN, 1, RW(p_xtmp+1), 1)
         call DAE2NLP(3, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT, RW(p_xtmp+1),
     1        iflag, IVAR, NNZA, RW(p_atmp+1), dummy, dummy, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-3tmp return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         p_rwend = p_xtmp
C
C     Compute Atmp * CSCALE / h
C
         call DAE2NLP(5, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT, XORIG,
     1        iflag, IVAR, NNZA, RW(p_atmp+1), CSCALE, VOUT, dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-5tmp return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT, XORIG,
     1        iflag, IVAR, NNZA, RW(p_atmp+1), CSCALE, VOUT(M+1), dummy,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-7tmp return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         p_rwend = p_atmp
C
c         p_xtmp  = p_rwend
c         p_rwend = p_xtmp + M
c         if( p_rwend.gt.LRW ) then
c            IERR = 98
c            goto 9999
c         endif
c         call DAE2NLP(5, NZ, NY, NU_OPT, NP_OPT, NCOL, NELE, TI,
C     1        ZINIT, XORIG,
c     1        iflag, IVAR, NNZA, RS(P_A+1), CSCALE, RW(p_xtmp+1),
c     2        dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
c     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
c     4        IERR)
c         if( IERR.gt.0 ) then
c            write(line,*) 'CONSTR: DAE2NLP-5 return IERR = ',IERR
c            call C_OUT(2,0,1,line)
c            goto 9999
c         endif
c         call DAE2NLP(7, NZ, NY, NU_OPT, NP_OPT, NCOL, NELE, TI,
C     1        ZINIT, XORIG,
c     1        iflag, IVAR, NNZA, RS(P_A+1), CSCALE, RW(p_xtmp+M+1),
c     2        dummy, LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
c     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
c     4        IERR)
c         if( IERR.gt.0 ) then
c            write(line,*) 'CONSTR: DAE2NLP-7 return IERR = ',IERR
c            call C_OUT(2,0,1,line)
c            goto 9999
c         endif
C
c         call DAXPY(NORIG, -1d0, RW(p_xtmp+1), 1, VOUT, 1)
         call DSCAL(NORIG, 1/h, VOUT, 1)

      else
CWEG
c         call DCOPY(NORIG, 0d0, 0, VOUT, 1)
c         goto 9999
c         do i = 1, M
c            write(*,*) 'LAM(',i,') = ',CSCALE(i)
c         enddo
         call DAE2NLP(11, NZ, NY, NU_OPT, NP_OPT, NCOL, NAC, NELE, TI,
     1        ZINIT, XORIG,
     1        iflag, IVAR, 0, dummy, VIN, VOUT, CSCALE,
     2        LRSD2N, RS(P_RSD2N+1), LISD2N, IS(P_ISD2N+1),
     3        LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     4        IERR)
         if( IERR.gt.0 ) then
            write(line,*) 'CONSTR: DAE2NLP-11 return IERR = ',IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif

C-------------------------------------------------------------------------------
C     End:   Compute Vout = (weighted Hessian of constraints) * Vin
C-------------------------------------------------------------------------------
      goto 9999

 9008 IERR = 8
      goto 9999

 9999 continue
      return
      end

