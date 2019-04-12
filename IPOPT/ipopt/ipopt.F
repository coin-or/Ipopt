C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      double precision function IPOPT(N, X, M, NLB, ILB, BNDS_L, NUB,
     1     IUB, BNDS_U, V_L, V_U, LAM, C, LRW, RW, LIW, IW, ITER, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3     DAT, IDAT, NARGS, ARGS, CARGS)
#ifdef BUILD_VSDLL
CDEC$ ATTRIBUTES DLLEXPORT::IPOPT
#endif
C
C*******************************************************************************
C
C    $Id: ipopt.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Main routine for IPOPT
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
CP   N         I    INT    total number of variables
CP   X        I/O   INT    I: user provided initial point
CP                         O: solution (hopefully)
CP                            (original order as in problem statement)
CP   M         I    INT    number of equality constraints
CP   NLB       I    INT    number of lower bounds
CP   ILB       I    INT    indices of lower bounds
CP                            (e.g. BNDS_L(i) is bound for X(ILB(i)) )
CP   BNDS_L    I    DP     values of lower bounds
CP   NUB       I    INT    number of upper bounds
CP   IUB       I    INT    indices of upper bounds
CP                            (e.g. BNDS_U(i) is bound for X(IUB(i)) )
CP   BNDS_U    I    DP     values of upper bounds
CP   V_L      I/O   DP     multipliers corresponding to the lower bounds
CP                            (ordered like BNDS_L) - set to zero for fixed
CP                            variables (input for warm start)
CP   V_U      I/O   DP     multipliers corresponding to the upper bounds
CP                            (ordered like BNDS_U) - set to zero for fixed
CP                            variables (input for warm start)
CP   LAM      I/O   DP     multipliers corresponding to the equality
CP                            constraints (input for warm start)
CP   C         O    DP     value of equality constraints at final point
CP   LRW       I    INT    length of RW
CP                         (if both LRW and LIW are zero, IPOPT allocates
CP                          memory on its own - if that option is compiled)
CP   RW        W    DP     DP work space
CP                            (internal usage: first RS, then RW)
CP   LIW       I    INT    length of IW
CP                         (if both LRW and LIW are zero, IPOPT allocates
CP                          memory on its own - if that option is compiled)
CP   IW        W    INT    INT work space
CP                            (internal usage: first IS, then IW)
CP   ITER      O    INT    number of iterations
CP   IERR      O    INT    Error code
CP                         =0: optimal solution found
CP                         =1: maximal number of iterations exceeded
CP                             (increase IMAXITER)
CP                         =2: step size becomes too small
CP                         =3: CPU time exceeded
CP                         =4: some algorithmic parameters wrong
CP                         =5: too few degrees of freedom
CP                         =6: it seems we have linear dependent
CP                             equality constraints
CP                         =7: some of the bounds are inconsistent
CP                         =8: error opening or reading a file
CP                         =9: Search direction getting too small
CP                        =10: Linear system becomes too ill-conditioned
CP                        =11: NaN or Inf occured
CP                        =13: Iterates becomes very large, seem to diverge
CP                        =15: Index in ILB or IUB out of range
CP                        =16: Point is (almost) feasible, but filter reverts
CP                             to restoration phase.
CP                        =17: Convergence to stationary point for
CP                             infeasibility
CP                             (Problem locally infeasible?)
CP                        =18: Restoration phase cannot further improve
CP                             feasibility
CP                        =20: The file 'STOP' has been detected in current
CP                             directory and caused IPOPT to abort
CP                        =96: Error during dynamic memory allocation
CP                        =97: An option was chosen for which a third party
CP                             component has not been compiled into excutable
CP                        =98: LRW too small
CP                        =99: LIW too small
CP                      other: Look at output... ;-)
CP   EV_F      I    EXT   Subroutine for objective function
CP   EV_C      I    EXT   Subroutine for constraints
CP   EV_G      I    EXT   Subroutine for gradient of objective fuction
CP   EV_A      I    EXT   Subroutine for Jacobian
CP   EV_H      I    EXT   Subroutine for Lagrangian Hessian
CP   EV_HLV    I    EXT   Subroutine for Lagrangian Hessian-vector products
CP   EV_HOV    I    EXT   Subroutine for objective Hessian-vector products
CP   EV_HCV    I    EXT   Subroutine for constraint Hessian-vector products
CP   DAT       P    DP    privat DP data for evaluation routines
CP                          will be passed untouched to your EVAL_* routines
CP   IDAT      P    INT   privat INT data for evaluation routines
CP                          will be passed untouched to your EVAL_* routines
CP   NARGS     I    INT    = 0 read the PARAMS.DAT file, if it can be
CP                             found, otherwise use default parameters
CP                         > 0 number of entries in ARGS, CARGS; if
CP                             the PARAMS.DAT file is present in current
CP                             directory, use the this in addition
CP                             (overwriting ARGS, CARGS)
CP                         < 0 -NARGS is number of entries in ARGS, CARGS
CP                             ignore the PARAMS.DAT file
CP   ARGS      I    DP     values to be assigned to corresponding
CP                         parameters listed in CARGS
CP   CARGS     I    C*(*)  names of parameters to be assigned
CP
CP   IPOPT     O    DP     return value: Objective function at final point
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
CS    ESTIMATE_WS
CS    IP_MALLOC
CS    IPOPTMAIN
CS    INITPARAMS
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
      integer N
      double precision X(N)
      integer M
      integer NLB
      integer ILB(NLB)
      double precision BNDS_L(NLB)
      integer NUB
      integer IUB(NUB)
      double precision BNDS_U(NUB)
      double precision V_L(NLB)
      double precision V_U(NUB)
      double precision LAM(M)
      double precision C(M)
      integer LRW
      double precision RW(*)
      integer LIW
      integer IW(*)
      integer ITER
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
      integer NARGS
      double precision ARGS(*)
      character*(*) CARGS(*)
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      character*80 line(4)
      double precision IPOPTMAIN
#ifdef USE_MALLOC
      integer IP_MALLOC, type, lrw1, liw1, i
#ifdef SIZEOF_INT_P_IS_4
      integer*4 ptr_rw, ptr_iw
#elif SIZEOF_INT_P_IS_8
      integer*8 ptr_rw, ptr_iw
#else
      Size of integer pointers unknown
#endif
#endif
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      IERR = 0
C
C     Call INITPARAMS to set the algorithmic parameters
C
      call INITPARAMS(IERR, NARGS, ARGS, CARGS)
      if( IERR.ne.0 ) then
         write(line,*)
     1        'ipopt: initparams ends with IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     If desired, do dynamic memory allocation
C
#ifdef USE_MALLOC
C
C     If required allocate memory
C
      ptr_rw = 0
      ptr_iw = 0
      call ESTIMATE_WS(N, X, M, NLB, NUB, 0.d0, lrw1, liw1, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1     DAT, IDAT)
      if( IERR.ne.0 ) then
         ITER = -1
         write(line,*) 'ipopt: ESTIMATE_WS returns IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( QPRINT.ge.2 ) then
         write(line,3) lrw1, liw1
 3       format(/,
     1'  Going to allocate double precision work space of size ',i16,/
     1'                             integer work space of size ',i16,/)
         call C_OUT(2,0,4,line)
      endif
      type = 1
      i = IP_MALLOC(type, lrw1, ptr_rw)
      if( i.ne.0 ) then
         write(line,*)
     1        'ipopt: Error allocating double precision work space: ',i
         IERR = 96
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      type = 2
      i = IP_MALLOC(type, liw1, ptr_iw)
      if( i.ne.0 ) then
         write(line,*)
     1        'ipopt: Error allocating integer work space: ',i
         IERR = 96
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      IPOPT = IPOPTMAIN(N, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U,
     1     V_L, V_U, LAM, C, lrw1, %val(ptr_rw), liw1, %VAL(ptr_iw),
     2     ITER, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV,
     2     EV_HCV, DAT, IDAT)

      type = 1
      i = IP_MALLOC(type, 0, ptr_rw)
      if( i.ne.0 ) then
         write(line,*)
     1        'ipopt: Error freeing double precision work space: ',i
         IERR = 96
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      type = 2
      i = IP_MALLOC(type, 0, ptr_iw)
      if( i.ne.0 ) then
         write(line,*)
     1        'ipopt: Error freeing integer work space: ',i
         IERR = 96
         call C_OUT(2,0,1,line)
         goto 9999
      endif
#else
      IPOPT=IPOPTMAIN(N, X, M, NLB, ILB, BNDS_L, NUB, IUB, BNDS_U, V_L,
     1     V_U, LAM, C, LRW, RW, LIW, IW, ITER, IERR, EV_F, EV_C, EV_G,
     2     EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
#endif
C
C     Close output file if it was opened
C
      if( QCNR.gt.0 ) call C_OUT(-2,0,0,' ')

 9999 continue
      return
      end

C*******************************************************************************
C
      double precision function IPOPTMAIN(NORIG, XORIG, M, NLBO, ILBO,
     1     BNDS_LO,NUBO, IUBO, BNDS_UO, V_LO, V_UO, LAM, C, LRW, RW,
     2     LIW, IW, ITER, IERR, EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV,
     3     EV_HOV, EV_HCV, DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: ipopt.F 574 2004-04-25 22:42:37Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Driver routine for IPOPT (without memory allocation)
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
CP   NORIG     I    INT    total number of variables (incl. fixed vars)
CP   XORIG    I/O   INT    I: user provided initial point (incl fixed vars)
CP                         O: solution (hopefully)
CP                            (original order as in problem statement)
CP   M         I    INT    number of equality constraints
CP   NLBO      I    INT    number of lower bounds (incl fixed vars)
CP   ILBO      I    INT    indices of lower bounds
CP                            (e.g. BNDS_LO(i) is bound for XORIG(ILBO(i)) )
CP   BNDS_LO   I    DP     values of lower bounds
CP   NUBO      I    INT    number of upper bounds (incl fixed vars)
CP   IUBO      I    INT    indices of upper bounds
CP                            (e.g. BNDS_UO(i) is bound for XORIG(IUBO(i)) )
CP   BNDS_UO   I    DP     values of upper bounds
CP   V_LO     I/O   DP     multipliers corresponding to the lower bounds
CP                            (ordered like BNDS_LO) - set to zero for fixed
CP                            variables (input for warm start)
CP   V_UO     I/O   DP     multipliers corresponding to the upper bounds
CP                            (ordered like BNDS_UO) - set to zero for fixed
CP                            variables (input for warm start)
CP   LAM      I/O   DP     multipliers corresponding to the equality
CP                            constraints (input for warm start)
CP   C         O    DP     value of equality constraints at final point
CP   LRW       I    INT    length of RW
CP   RW        W    DP     DP work space
CP                            (internal usage: first RS, then RW)
CP   LIW       I    INT    length of IW
CP   IW        W    INT    INT work space
CP                            (internal usage: first IS, then IW)
CP   ITER      O    INT    number of iterations
CP   IERR      O    INT    Exit code
CP   EV_F      I    EXT    Subroutine for objective function
CP   EV_C      I    EXT    Subroutine for constraints
CP   EV_G      I    EXT    Subroutine for gradient of objective fuction
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
CS    GET_IFIX
CS    INIT_MEM
CS    INIT_BND
CS    MAINLOOP
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
      include 'TIMER.INC'
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer NORIG
      double precision XORIG(NORIG)
      integer M
      integer NLBO
      integer ILBO(NLBO)
      double precision BNDS_LO(NLBO)
      integer NUBO
      integer IUBO(NUBO)
      double precision BNDS_UO(NUBO)
      double precision V_LO(NLBO)
      double precision V_UO(NUBO)
      double precision LAM(M)
      double precision C(M)
      integer LRW
      double precision RW(LRW)
      integer LIW
      integer IW(LIW)
      integer ITER
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
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer KCONSTR(6)
      save    KCONSTR

      integer LRS_END, LIS_END, P_IFIX, P_X, P_BNDSL, P_BNDSU
      save    LRS_END, LIS_END, P_IFIX, P_X, P_BNDSL, P_BNDSU
      integer P_ILB, P_IUB, P_IVAR, P_CSCALE
      save    P_ILB, P_IUB, P_IVAR, P_CSCALE
      integer N, NIND, NLB, NUB, NFIX, NZORIG
      save    N, NIND, NLB, NUB, NFIX, NZORIG

      double precision stime, ftime, f
      integer p_rwend, p_iwend, p_vcopy, p_ibo1, p_ilbo1, p_iubo1
      integer i, idummy

      character*80 line(6)
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      call TIMER(stime)
      IERR = 0
      f = 0.d0
C
C     Initialize work space pointers
C
      p_rwend = 0
      p_iwend = 0
      LRS_END = 0
      LIS_END = 0
C
C     initialize iteration counter
C
      ITER = 0
C
C     Find fixed variables
C
      P_IFIX  = LIS_END
      LIS_END = P_IFIX  + NORIG
      p_iwend = LIS_END

      p_ilbo1 = p_iwend
      p_iubo1 = p_ilbo1 + NORIG
      p_iwend = p_iubo1 + NORIG
      if( p_iwend.gt.liw ) then
         IERR = 99
         goto 9999
      endif
      call GET_IFIX(NORIG, XORIG, NLBO, ILBO, BNDS_LO,
     1              NUBO, IUBO, BNDS_UO, NFIX, IW(P_IFIX+1),
     2              IW(p_ilbo1+1),IW(p_iubo1+1), IERR)
      if( IERR.ne.0 ) then
         if( IERR.eq.7 ) goto 8000
         write(line,*) 'ipopt: GET_IFIX returns IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      LIS_END = P_IFIX + NFIX
      p_iwend = LIS_END
C
C     Initialize some sizes
C
      N    = NORIG - NFIX
      NIND = N - M
      NLB  = NLBO - NFIX
      NUB  = NUBO - NFIX
C
C     Print problem size
C
      write(line,710) NORIG, NFIX, M, NLB, NUB
 710  format(/,'Number of variables           : ',i8,
     1       /,'   of which are fixed         : ',i8,
     1       /,'Number of constraints         : ',i8,
     1       /,'Number of lower bounds        : ',i8,
     1       /,'Number of upper bounds        : ',i8)
      call C_OUT(2,0,6,line)

      if( NIND.lt.0 ) then
         write(line,*)
     1        'Error: Number of degrees of freedom is NIND = ',NIND
         call C_OUT(2,0,1,line)
         IERR = 5
         goto 8000
      endif
C
C     Determine demand of storage and work space
C
      INMEMCHECK = 0
      call INIT_MEM(NORIG, XORIG, N, NIND, M, NLB, NUB, NZORIG, idummy,
     1     KCONSTR, LRS_END, LIS_END, LRW-p_rwend, RW(p_rwend+1),
     1     LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'ipopt: INIT_MEM returns IERR = ', IERR
         call C_OUT(2,0,1,line)
        goto 9999
      endif
      if( LRS_END.gt.LRW ) then
         IERR = 98
         call C_OUT(2,0,1,'Not enough double precision storage space:')
         write(line,*) 'LRS_END = ',LRS_END,', but LRW = ',lrw
         call C_OUT(2,0,1,line)
         goto 8000
      elseif( LIS_END.gt.LIW ) then
         IERR = 99
         call C_OUT(2,0,1,'Not enough integer storage space:')
         write(line,*) 'LIS_END = ',LIS_END,', but LIW = ',liw
         call C_OUT(2,0,1,line)
         goto 8000
      endif
C
C     Now reserve storage space for this subroutine
C
      P_X     = LRS_END
      P_BNDSL = P_X     + N
      P_BNDSU = P_BNDSL + NLB
      LRS_END = P_BNDSU + NUB
      P_ILB   = LIS_END
      P_IUB   = P_ILB   + NLB
      P_IVAR  = P_IUB   + NUB
      LIS_END = P_IVAR  + N
      if( QSCALE.ge.3 ) then
         P_CSCALE  = LRS_END
         LRS_END   = P_CSCALE + M + N
      elseif( QSCALE.ge.2 ) then
         P_CSCALE  = LRS_END
         LRS_END   = P_CSCALE + M
      else
         P_CSCALE  = LRS_END
         LRS_END   = P_CSCALE + 1
      endif
C
C     The rest in RW and IW is now available as work space
C
      p_rwend = LRS_END
      p_iwend = LIS_END
      if( MEMDBG ) then
         write(line,1)'ipopt', LRS_END, LIS_END
 1       format('MEMDBG - ',a20,': LRS_END = ',i12,' LIS_END = ',i12)
         call C_OUT(1,0,1,line)
         write(line,2)'ipopt', LRW, LIW
 2       format('MEMDBG - ',a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 8000
      elseif( p_iwend.gt.LIW ) then
         IERR = 99
         goto 8000
      endif
C
C     If required, compute scaling factors
C
      if( QSCALE.ne.0 ) then
         call GET_SCALE(ITER, NORIG, XORIG, N, RW(P_X+1), NIND, M,
     1        IW(P_IVAR+1), NFIX, IW(P_IFIX+1), NLBO, ILBO, BNDS_LO,
     1        NUBO, IUBO, BNDS_UO, RW(P_CSCALE+1), KCONSTR, LRS_END,
     1        RW, LIS_END, IW, LRW-p_rwend, RW(p_rwend+1),
     1        LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'ipopt: get_scale returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 8000
         endif
      else
         RW(P_CSCALE+1) = 1.d0
      endif
C
C     Now initialize X, BNDS_L, BNDS_U, and rescale VO_L, VO_U, LAM if necessary
C
      call INIT_BND(NORIG, N, NIND, M, NFIX, IW(P_IFIX+1), XORIG,
     1     RW(P_X+1), NLBO, ILBO, BNDS_LO, V_LO, NLB, IW(P_ILB+1),
     1     RW(P_BNDSL+1), NUBO, IUBO, BNDS_UO, V_UO, NUB, IW(P_IUB+1),
     1     RW(P_BNDSU+1), IW(P_IVAR+1), RW(P_CSCALE+1), LAM,
     1     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1), IERR)
      if( IERR.ne.0 ) then
         write(line,*) 'ipopt: init_bnd returns IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 8000
      endif
C
C     After all this initialization call algorithm
C
      call MAINLOOP(ITER, N, NIND, M, NORIG, XORIG,
     1     IW(P_IVAR+1), RW(P_X+1), NFIX, IW(P_IFIX+1),
     1     NLB, IW(P_ILB+1), NUB, IW(P_IUB+1),
     1     RW(P_BNDSL+1), RW(P_BNDSU+1), NZORIG,
     1     V_LO, V_UO, LAM, f, C, RW(P_CSCALE+1), KCONSTR,
     1     LRS_END, LRS_END, RW, LIS_END, LIS_END, IW,
     1     LRW-p_rwend, RW(p_rwend+1), LIW-p_iwend, IW(p_iwend+1),
     1     IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C     Change IERR flag in case we actually solved a feasibility problem
C
      if( NIND.eq.0 .and. IERR.eq.16 ) IERR = 0
C
C     Correct the position of the bound multipliers in V_LO and V_UO
C
      p_vcopy = p_rwend
      p_rwend = p_vcopy + max(NLB,NUB)
      if( p_rwend.gt.LRW ) then
         IERR = 98
         goto 9999
      endif
      p_ibo1  = p_iwend
      p_iwend = p_ibo1 + NORIG
      if( p_iwend.gt.LIW ) then
         IERR = 99
         goto 9999
      endif
      call DCOPY(NLB, V_LO, 1, RW(p_vcopy+1), 1)
      if( NFIX.gt.0 ) then
         call DCOPY(NLBO, 0.d0, 0, V_LO, 1)
      endif
      do i = 1, NLBO
         IW(p_ibo1+ILBO(i)) = i
      enddo
      do i = 1, NLB
         V_LO(IW(p_ibo1+IW(P_IVAR+IW(P_ILB+i)))) = RW(p_vcopy+i)
      enddo
      call DCOPY(NUB, V_UO, 1, RW(p_vcopy+1), 1)
      if( NFIX.gt.0 ) then
         call DCOPY(NUBO, 0.d0, 0, V_UO, 1)
      endif
      do i = 1, NUBO
         IW(p_ibo1+IUBO(i)) = i
      enddo
      do i = 1, NUB
         V_UO(IW(p_ibo1+IW(P_IVAR+IW(P_IUB+i)))) = RW(p_vcopy+i)
      enddo

 8000 continue
C
C     Write status line
C
      line(1) = ' '
      if( IERR.eq.0 ) then
         line(2) = 'EXIT: OPTIMAL SOLUTION FOUND'
      elseif( IERR.eq.1 ) then
         line(2) = 'EXIT: Maximal number of iterations exceeded'
      elseif( IERR.eq.2 ) then
         line(2) = 'EXIT: Step size becomes too small'
      elseif( IERR.eq.3 ) then
         line(2) = 'EXIT: Allowed CPU time exceeded'
      elseif( IERR.eq.4 ) then
         line(2) =
     1  'EXIT: Some algorithmic parameters are wrong or inconsistent'
      elseif( IERR.eq.5 ) then
         line(2) = 'EXIT: Too few degrees of freedom'
      elseif( IERR.eq.6 ) then
         line(2) =
     1      'EXIT: Equality constraints seem to be linearly dependent'
      elseif( IERR.eq.7 ) then
         line(2) = 'EXIT: Some of the bounds are inconsistent'
      elseif( IERR.eq.8 ) then
         line(2) = 'EXIT: Error opening or reading some file'
      elseif( IERR.eq.9 ) then
         line(2) =
     1  'EXIT: Search direction becomes too small'
      elseif( IERR.eq.10 ) then
         line(2) =
     1        'EXIT: Linear system becomes too ill-conditioned'
      elseif( IERR.eq.11 ) then
         line(2) =
     1        'EXIT: NaN or Inf occured'
      elseif( IERR.eq.13 ) then
         line(2) =
     1        'EXIT: Iterates become very large (diverging?)'
      elseif( IERR.eq.15 ) then
         line(2) = 'EXIT: Index in ILB or IUB out of range'
      elseif( IERR.eq.16 ) then
         line(2) =
     1 'EXIT: Point is (almost) feasible, but restoration phase called'
      elseif( IERR.eq.17 ) then
         line(2) =
     1        'EXIT: Convergence to stationary point for infeasibility'
      elseif( IERR.eq.18 ) then
         line(2) =
     1      'EXIT: Restoration phase cannot further improve feasibility'
      elseif( IERR.eq.20 ) then
         line(2) = 'EXIT: STOP file detected'
      elseif( IERR.eq.96 ) then
         line(2) = 'EXIT: Error during dynamic memory allocation'
      elseif( IERR.eq.97 ) then
         line(2) =
     1        'EXIT: Third party component missing for chosen options'
      elseif( IERR.eq.98 ) then
         line(2) = 'EXIT: Size of double precision work space too small'
      elseif( IERR.eq.99 ) then
         line(2) = 'EXIT: Size of integer work space too small'
      else
         write(line, 1000) IERR
 1000    format(/,'EXIT: IERR =',i5,'; check output.')
      endif
      call C_OUT(2,0,2,line)


      if( QPRINT.ge.0 ) then
         call C_OUT(1,0,1,' ')
         write(line,*) 'TIME_GET_STEP_FULL = ', TIME_GET_STEP_FULL
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_YPY           = ', TIME_YPY
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_BB            = ', TIME_BB
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_CG            = ', TIME_CG
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_HV            = ', TIME_HV
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_EXACTW        = ', TIME_EXACTW
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_ZWZY_BACKS    = ', TIME_ZWZY_BACKS
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_ZWZY_EGVAL    = ', TIME_ZWZY_EVALA
         call C_OUT(1,0,1,line)
         write(line,*) 'TIME_PZ_CHOL       = ', TIME_PZ_CHOL
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_CG           = ', COUNT_CG
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_HV           = ', COUNT_HV
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_RESTO_CALL   = ', COUNT_RESTO_CALL
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_RESTO_ITER   = ', COUNT_RESTO_ITER
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_TRON_CG      = ', COUNT_TRON_CG
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_NEG_CURV     = ', COUNT_NEG_CURV
         call C_OUT(1,0,1,line)
         write(line,*) 'COUNT_DEPCON       = ', COUNT_DEPCON
         call C_OUT(1,0,1,line)
         call C_OUT(1,0,1,' ')
      endif

      if( PRINTTIMING ) then
         write(*,1001) TIME_GET_STEP_FULL
 1001 format(/,'CPU seconds spent in GET_STEP_FULL                  = ',
     1        f15.4)
      endif

 9999 continue
      IPOPTMAIN = f
      call TIMER(ftime)
      write(line,9998) ftime-stime
 9998 format(/,'CPU seconds spent in IPOPT and function evaluations = ',
     1     f15.4,/)
      call C_OUT(2,0,3,line)
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine IPOPT_WS(N, X, M, NLB, NUB, NZA, NZH, LRW, LIW,
     1     FILLINFACT, EV_H, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, NZH, LRW, LIW
      external EV_H
      double precision DAT(*)
      integer IDAT(*)
      double precision FILLINFACT, X(N)
      integer lrw1, liw1
      character*80 line

      call INIT_BND_WS(N, M, NLB, NUB, NZA, LRW, LIW)
      if( QSCALE.ne.0 ) then
         call GET_SCALE_WS(N, M, NLB, NUB, NZA, NZH, lrw1, liw1, DAT,
     1        IDAT)
         LRW = max(LRW,lrw1)
         LIW = max(LIW,liw1)
      endif
      call MAINLOOP_WS(N, X, M, NLB, NUB, NZA, lrw1, liw1,
     1     FILLINFACT, EV_H, DAT, IDAT)
      LRW = max(LRW, lrw1, NLB, NUB)
      LIW = max(LIW, liw1, N)

      if( QPRINT.ge.4 ) then
         write(line,1000) 'ipopt_ws',LRW,LIW
 1000    format(a20,': LRW = ',i12,' LIW = ',i12)
         call C_OUT(1,0,1,line)
      endif

      return
      end

