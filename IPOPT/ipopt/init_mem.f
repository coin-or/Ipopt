C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine INIT_MEM(NORIG, XORIG, N, NIND, M, NLB, NUB, NZORIG,
     1     NNZHORIG, KCONSTR, LRS_END, LIS_END, LRW, RW, LIW, IW, IERR,
     2     EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     3     DAT, IDAT)

C
C*******************************************************************************
C
C    $Id: init_mem.f 655 2004-10-05 17:23:15Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Main (outer) loop of algorithm
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
CP   XORIG     I    DP     initial values of variables (incl. fixed vars)
CP   N         I    INT    number of variables (without fixed)
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of constraints
CP   NLB       I    INT    number of lower bounds (excluding fixed vars)
CP   NUB       I    INT    number of upper bounds (excluding fixed vars)
CP   NZORIG    I    INT    number of nonzeros in Jacobian of constraints
CP                            (including rows to fixed variables!)
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP   LRS_END  I/O   INT    last used reserved entry in RS
CP   LRS_END  I/O   INT    last used reserved entry in IS
CP   LRW       I    INT    length of RW
CP   RW       I/O   DP     can be used as DP work space but content will be
CP                            changed between calls
CP   LIW       I    INT    length of IW
CP   IW       I/O   INT    can be used as INT work space but content will be
CP                            changed between calls
CP   IERR      O    INT    =0: everything OK
CP                         >0: Error occured; abort optimization
CP                         <0: Warning; message to user
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
CS    CONSTR
CS    MAINLOOP
CS    LINESEARCH
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
      integer N
      integer NIND
      integer M
      integer NLB
      integer NUB
      integer NZORIG
      integer NNZHORIG
      integer KCONSTR(6)
      integer LRS_END
      integer LIS_END
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
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer lrs_constr, lis_constr
      integer idummy
      double precision dummy
      logical ldummy
      character*1 cdummy
      character*80 line

      integer LRW_CONSTR, LIW_CONSTR
      common /CONSTRWS/ LRW_CONSTR, LIW_CONSTR
      save   /CONSTRWS/
C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C

C
C     Initialize CPU times
C
      TIME_BB            = 0.d0
      TIME_CG            = 0.d0
      TIME_YPY           = 0.d0
      TIME_EXACTW        = 0.d0
      TIME_ZWZY_BACKS    = 0.d0
      TIME_ZWZY_EVALA    = 0.d0
      TIME_PZ_CHOL       = 0.d0
      TIME_HV            = 0.d0
      TIME_GET_STEP_FULL = 0.d0
      COUNT_CG           = 0
      COUNT_RESTO_ITER   = 0
      COUNT_NEG_CURV     = 0
      COUNT_RESTO_CALL   = 0
      COUNT_TRON_CG      = 0
      COUNT_HV           = 0
      COUNT_DEPCON       = 0
C
C     call subroutines to initialize their pointers for the storage space
C

C
C     CONSTR
C
      if( M.ne.0 ) then
         lrw_constr = LRW
         liw_constr = LIW
         call CONSTR(0, idummy, N, NIND, M, idummy, idummy, idummy,
     1        NORIG, XORIG, dummy, dummy, dummy, NZORIG,
     2        idummy, lrs_constr, dummy, lis_constr, idummy,
     3        LRW_CONSTR, RW, LIW_CONSTR, IW, IERR, EV_F, EV_C, EV_G,
     5        EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: CONSTR returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         KCONSTR(1) = lrs_constr
         KCONSTR(3) = lis_constr
         KCONSTR(5) = LRW_CONSTR
         KCONSTR(6) = LIW_CONSTR

         KCONSTR(2) = LRS_END
         LRS_END    = KCONSTR(2) + lrs_constr
         KCONSTR(4) = LIS_END
         LIS_END    = KCONSTR(4) + lis_constr
      else
         NZORIG = 0
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'init_mem: After CONSTR LRS_END = ',LRS_END,
     1        ' LIS_END = ',LIS_END
         call C_OUT(1,0,1,line)
      endif
      if( INMEMCHECK.eq.0 ) then
         write(line,1000) NZORIG
 1000    format('Number of nonzeros in Jacobian: ',i8)
         call C_OUT(2,0,1,line)
      endif
C
C     MAINLOOP
C
      call MAINLOOP(-1, N, NIND, M, NORIG, dummy,
     1     idummy, dummy, idummy, idummy, NLB, idummy, NUB,
     2     idummy, dummy, dummy, NZORIG, dummy, dummy,
     1     dummy, dummy, dummy, dummy, idummy, idummy,
     3     LRS_END, dummy, idummy, LIS_END, idummy,
     4     idummy, dummy, idummy, idummy, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'init_mem: MAINLOOP returns IERR = ', IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'init_mem: After MAINLOOP LRS_END = ',LRS_END,
     1        ' LIS_END = ',LIS_END
         call C_OUT(1,0,1,line)
      endif
C
C     LINESEARCH
C
      if( abs(QMERIT).eq.1 .or. abs(QMERIT).eq.2 .or. QMERIT.eq.0 ) then
         call LINESEARCH(-1, N, NIND, M, dummy, idummy, NLB, idummy,
     1        NUB, idummy, dummy, dummy, dummy, dummy, dummy,
     2        dummy, dummy, dummy, dummy, dummy, dummy, NORIG,
     3        dummy, dummy, dummy, dummy, dummy,
     4        dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     2        dummy, dummy, dummy, ldummy, dummy, dummy, idummy,
     5        cdummy, dummy, ldummy, idummy, idummy,
     5        idummy, LRS_END, dummy, idummy, LIS_END,
     6        idummy, idummy, dummy, idummy, idummy, IERR, EV_F, EV_C,
     7        EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: LINESEARCH returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( QPRINT.ge.2 ) then
            write(line,*) 'init_mem: After LINESEARCH LRS_END = ',
     1           LRS_END,' LIS_END = ',LIS_END
            call C_OUT(1,0,1,line)
         endif
      endif
      if( abs(QMERIT).ge.4 ) then
         call FILTER(-1, N, NIND, M, dummy, idummy, idummy, idummy,
     1        NLB, idummy, NUB, idummy, dummy, dummy, dummy, dummy,
     2        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     3        NORIG, dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     4        dummy, dummy, dummy, dummy, dummy, idummy, idummy,
     5        dummy, cdummy, cdummy, idummy, idummy, ldummy, dummy,
     5        dummy, dummy, dummy, dummy, dummy, dummy, idummy, LRS_END,
     6        dummy, LIS_END, idummy, idummy, dummy, idummy, idummy,
     7        IERR, EV_F, EV_C, EV_G, EV_A,
     5        EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: FILTER returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( QPRINT.ge.2 ) then
            write(line,*) 'init_mem: After FILTER LRS_END = ',
     1           LRS_END,' LIS_END = ',LIS_END
            call C_OUT(1,0,1,line)
         endif
      endif
      if( abs(QMERIT).ge.4 ) then
         if( abs(QRESTO).eq.1 ) then
            call RESTO_FILTER(-1, N, NIND, M, NORIG, XORIG,
     1           dummy, idummy, idummy, idummy, NLB, idummy, NUB,
     1           idummy, dummy, dummy, dummy, dummy, dummy, dummy,
     1           idummy, dummy, dummy, dummy, dummy, dummy, dummy,
     1           dummy, dummy, dummy, dummy, dummy,
     1           dummy, dummy, dummy, dummy, dummy,
     1           dummy, dummy, dummy, dummy, dummy, dummy,
     1           idummy, cdummy, dummy, dummy, idummy, dummy,
     1           dummy, dummy, dummy, dummy, dummy,
     1           idummy, LRS_END, dummy, LIS_END, idummy,
     1           idummy, dummy, idummy, idummy, IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         endif
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: RESTO_FILTER returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( QPRINT.ge.2 ) then
            write(line,*) 'init_mem: After RESTO_FILTER LRS_END = ',
     1           LRS_END,' LIS_END = ',LIS_END
            call C_OUT(1,0,1,line)
         endif
         if( QRESTO.eq.2 ) then
            call RESTO_TRON(-1, N, NIND, M, NORIG, XORIG,
     1           dummy, idummy, idummy, idummy, NLB, idummy, NUB,
     1           idummy, dummy, dummy, dummy, dummy, dummy, NZORIG,
     1           dummy, dummy, dummy, dummy, dummy,
     1           dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     1           dummy, dummy, dummy,
     1           idummy, LRS_END, dummy, LIS_END, idummy,
     1           idummy, dummy, idummy, idummy, IERR, EV_F, EV_C, EV_G,
     5           EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
         endif
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: RESTO_TRON returns IERR = ', IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( QPRINT.ge.2 ) then
            write(line,*) 'init_mem: After RESTO_TRON LRS_END = ',
     1           LRS_END,' LIS_END = ',LIS_END
            call C_OUT(1,0,1,line)
         endif
      endif
C
C     GET_STEP_FULL
C
      if( QFULL.ne.0 ) then
         call GET_STEP_FULL(-1, NORIG, N, NIND, M, NZORIG, dummy,
     1        XORIG, dummy, dummy, NLB, idummy, NUB, idummy, NNZHORIG,
     1        idummy, idummy, dummy, dummy,
     1        dummy, dummy, dummy, dummy, dummy, dummy,
     1        dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     1        dummy, dummy, dummy, dummy, cdummy, ldummy, idummy, dummy,
     1        dummy, dummy, dummy, dummy,
     1        dummy, idummy, idummy, dummy, dummy, idummy, idummy,
     1        idummy, ldummy, idummy, idummy, LRS_END, dummy, idummy,
     1        LIS_END, idummy, idummy, dummy, idummy, idummy, IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1        DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: GET_STEP_FULL returns IERR = ',
     1           IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
         if( QPRINT.ge.2 ) then
            write(line,*) 'init_mem: After GET_STEP_FULL LRS_END = ',
     1           LRS_END,' LIS_END = ',LIS_END
            call C_OUT(1,0,1,line)
         endif
      else
         NNZHORIG = 0
      endif
C
C     GET_SCALE (just to store NNZA and NNZH)
C
      if( QSCALE.ge.3 ) then
         call GET_SCALE(-1, NORIG, dummy, N, dummy, NZORIG, M, idummy,
     1        nnzhorig, idummy, idummy, idummy, dummy, idummy, idummy,
     1        dummy, dummy,
     1        idummy, idummy, dummy, idummy, idummy, idummy, dummy,
     1        idummy, idummy, IERR,
     2        EV_F, EV_C, EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV,
     1        DAT, IDAT)
         if( IERR.ne.0 ) then
            write(line,*) 'init_mem: GET_SCALE returns IERR = ',
     1           IERR
            call C_OUT(2,0,1,line)
            goto 9999
         endif
      endif
C
C     ERROR (just for initialization of LASTITER)
C
      call OPTERROR(N, NIND, M, dummy, dummy, dummy, dummy, dummy,
     1     dummy, NLB, idummy, NUB, idummy, dummy, dummy,
     1     dummy, dummy, dummy, dummy, dummy, dummy, dummy, cdummy,
     2     dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy,
     1     dummy, dummy, -1, idummy, idummy, idummy, NORIG, dummy,
     1     dummy, KCONSTR, idummy, dummy, idummy, idummy,
     1     idummy, dummy, idummy, idummy, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'init_mem: OPTERROR returns IERR = ',
     1        IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'init_mem: After OPTERROR LRS_END = ',
     1        LRS_END,' LIS_END = ',LIS_END
         call C_OUT(1,0,1,line)
      endif
C
C     GET_HV
C
      call GET_HV(-1, N, NIND, idummy, idummy, dummy, idummy, NORIG,
     1     XORIG, NLB, idummy, NUB, idummy, dummy, dummy,
     2     dummy, M, dummy, dummy, dummy,
     2     KCONSTR, idummy, dummy, idummy, idummy, idummy, dummy,
     3     idummy, idummy, IERR, EV_F, EV_C, EV_G, EV_A,
     5     EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.ne.0 ) then
         write(line,*) 'init_mem: GET_HV returns IERR = ',
     1        IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
      if( QPRINT.ge.2 ) then
         write(line,*) 'init_mem: After GET_HV = ',
     1        LRS_END,' LIS_END = ',LIS_END
         call C_OUT(1,0,1,line)
      endif
C
 9999 continue
      return
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine CONSTR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)
      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)
      integer LRW_CONSTR, LIW_CONSTR
      common /CONSTRWS/ LRW_CONSTR, LIW_CONSTR
      save   /CONSTRWS/

      LIW = LIW_CONSTR
      LRW = LRW_CONSTR

      return
      end

