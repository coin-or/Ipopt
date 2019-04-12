C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine PARTITION(TASK, ITER, N, NIND, M, NFIX, IFIX, IVAR,
     1     NORIG, XORIG, CSCALE, KCONSTR, LRS, RS,
     2     LIS, IS, LRW, RW, LIW, IW, IERR, EV_F, EV_C,
     3     EV_G, EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
C
C*******************************************************************************
C
C    $Id: partition.f 531 2004-03-11 01:31:07Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Find out partition of variables into dependent and independent
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
CP   TASK      I    INT    =0: Let MA28 choose the basis
CP                         =1: Choose first variables as independent
CP                         =2: Choose last variables as independent
CP                         =3: Read independent variables from file
CP                         =4: Read dependent variables from file
CP   ITER      I    INT    iteration counter (if 0, this is first call)
CP   N         I    INT    number of (free) variables
CP   NIND      I    INT    number of independent variables
CP   M         I    INT    number of equality constraints / dependent variables
CP   NFIX      I    INT    number of fixed variables
CP   IFIX      I    INT    specifies variables that are fixed by bounds:
CP                            i = 1..NORIG-N   XORIG(IFIX(i)) is fixed
CP                            (assumed to be in increasing order)
CP   IVAR      O    INT    information about partitioning
CP                            i = 1..M      XORIG(IVAR(i)) dependent
CP                            i = (M+1)..N  XORIG(IVAR(i)) independent
CP                            Note: fixed variables do not occur in IVAR
CP   NORIG     I    INT    number of all variables including fixed vars
CP   XORIG     I    DP     (only for TASK = 0)
CP                            actual iteration (including fixed entries)
CP   CSCALE    I    DP     scaling factors for constraints
CP   KCONSTR   I    INT    KCONSTR(1): LRS for CONSTR
CP                         KCONSTR(2): P_LRS for CONSTR
CP                         KCONSTR(3): LIS for CONSTR
CP                         KCONSTR(4): P_LIS for CONSTR
CP                         KCONSTR(5): LRW for CONSTR
CP                         KCONSTR(6): LIW for CONSTR
CP   LRS       I    INT    total length of RS
CP   RS       I/O   DP     DP storage space (all!)
CP   LIS       I    INT    total length of IS
CP   IS       I/O   INT    INT storage space (all!)
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
      integer TASK
      integer ITER
      integer N
      integer NIND
      integer M
      integer NFIX
      integer IFIX(NFIX)
      integer IVAR(N)
      integer NORIG
      double precision XORIG(NORIG)
      double precision CSCALE(*)
      integer KCONSTR(6)
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
      integer p_iwend, p_rwend, p_xcor
      integer i, j, lind, ldep, lfix, idummy, icor, lx
      double precision dummy
      logical fixed

      character*80 line(6)
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
C
C     If no equality constraints, all variables are independent
C
      if( M.eq.0 ) goto 200
C
C     If there are no degrees of freedom, all variables are dependent
C
      if( NIND.eq.0 ) goto 300
C
C     Sanity check
C
      if( QSCALE.ge.3 ) then
         call C_OUT(2,0,1,'partition: Basis selection not properly imple
     1mented for this scaling option.')
         IERR = 4
         goto 9999
      endif
C
C     Depending on TASK select the independent variables
C
      goto (100, 200, 300, 400, 500) TASK + 1
      call C_OUT(2,0,1,'partition: Invalid flag for basis selection.')
      IERR = 4
      goto 9999

 100  continue
C
C     Start: Let MA28 choose the basis   (THIS IS IN CONSTR!!!!)
C
      call CONSTR(1, ITER, N, NIND, M, IVAR, NFIX, IFIX,
     1            NORIG, XORIG, CSCALE, dummy, dummy, idummy, idummy,
     2            KCONSTR(1), RS(KCONSTR(2)+1), KCONSTR(3),
     4            IS(KCONSTR(4)+1), LRW-p_rwend, RW(p_rwend+1),
     5            LIW-p_iwend, IW(p_iwend+1), IERR, EV_F, EV_C, EV_G,
     5            EV_A, EV_H, EV_HLV, EV_HOV, EV_HCV, DAT, IDAT)
      if( IERR.lt.0 ) then
         write(line,*) 'partition: Warning in CONSTR-1, IERR = ',IERR
         call C_OUT(2,0,1,line)
      elseif( IERR.ne.0 ) then
         write(line,*) 'partition: Error in CONSTR-1, IERR = ',IERR
         call C_OUT(2,0,1,line)
         goto 9999
      endif
C
C     End:   Let MA28 choose the basis
C
      goto 1000
C
C     Start: Choose first variables as independent
C
 200  continue
      do i = 1, M
         IVAR(i) = NIND + i
      enddo
      do i = 1, NIND
         IVAR(M+i) = i
      enddo
C
C     End:   Choose first variables as independent
C
      goto 1000
C
C     Start: Choose last variables as independent
C
 300  continue
      do i = 1, N
         IVAR(i) = i
      enddo
C
C     End:   Choose last variables as independent
C
      goto 1000
C
C     Start: Read independent variables from file 'BASIS.DAT'
C
C            It is assumed that they are ordered increasingly!
C            Numbering is according to XORIG and no fixed
C            variables are allowed to be mentioned!
C
 400  continue
      open(80,file='BASIS.DAT',status='old',err=9008)
      do i = 1, NIND
         read(80,410,end=9009,err=9009) IVAR(M+i)
 410     format(i15)
      enddo
      lfix = 1
      lind = M+1
      ldep = 1
      do i = 1, NORIG
         fixed = .false.
         if( lfix.le.NFIX ) then
            if( IFIX(lfix).eq.i ) then
               fixed = .true.
               lfix = lfix + 1
            endif
         endif
         if( .not.fixed ) then
            if( lind.le.N ) then
               if( IVAR(lind).eq.i ) then
                  lind = lind + 1
                  goto 420
               endif
            endif
            IVAR(ldep) = i
            ldep = ldep + 1
         endif
 420     continue
      enddo
C
C     End:   Read independent variables from file
C
      goto 1000
C
C     Start: Read dependent variables from file 'BASIS.DAT'
C
C            It is assumed that they are ordered increasingly!
C            Numbering is according to XORIG and no fixed
C            variables are allowed to be mentioned!
C
 500  continue
      open(80,file='BASIS.DAT',status='old',err=9008)
      do i = 1, M
         read(80,510,end=9009,err=9009) IVAR(i)
 510     format(i15)
      enddo
      lfix = 1
      lind = M+1
      ldep = 1
      do i = 1, NORIG
         fixed = .false.
         if( lfix.le.NFIX ) then
            if( IFIX(lfix).eq.i ) then
               fixed = .true.
               lfix = lfix + 1
            endif
         endif
         if( .not.fixed ) then
            if( ldep.le.M ) then
               if( IVAR(ldep).eq.i ) then
                  ldep = ldep + 1
                  goto 520
               endif
            endif
            IVAR(lind) = i
            lind = lind + 1
         endif
 520     continue
      enddo
C
C     End:   Read dependent variables from file
C
      goto 1000
C
C     Need to take care of fixed variables
C     (this assumes IFIX to be in increasing order)
C
 1000 continue
      if( NFIX.ne.0 ) then
C
C     Determine for each entry in partition, how much the corresponding entry
C     in IVAR has to be corrected
C
         p_xcor  = p_iwend
         p_iwend = p_xcor + N
         if( p_iwend.gt.LIW ) then
            IERR = 99
            goto 9999
         endif
         lfix = 1
         icor = 0
         lx = 0
         do i = 1, NORIG
            if( i.eq.IFIX(lfix) ) then   ! new fixed element
               icor = icor + 1
               lfix = lfix + 1
               if( lfix.gt.NFIX ) then
                  do j = 1, NORIG-i
                     IW(p_xcor+j+lx) = icor
                  enddo
                  goto 1010
               endif
            else
               lx = lx + 1
               IW(p_xcor+lx) = icor
            endif
         enddo
 1010    continue
C
C     Now correct IVAR
C
         do i = 1, N
            IVAR(i) = IVAR(i) + IW(p_xcor+IVAR(i))
         enddo
C
         p_iwend = p_xcor
C
      endif
C
C     Write partition to file
C
 1100 continue
      if( QCNR.gt.0 .and. QPRINT.ge.2 ) then
         write(line,700) ITER
 700     format(/,'   Information about partition in ITER ',i5,':',//,
     1          '   Dependent variables:',/)
         call C_OUT(1,2,5,line)
         do i = 1, M
            write(line,730) IVAR(i), XORIG(IVAR(i)), i
 730        format(' XORIG(',i5,') = ',d15.8,
     1             '  is   dependent number ',i5)
            call C_OUT(1,2,1,line)
         enddo
         write(line,720)
 720     format(/,'   Independent variables:',/)
         call C_OUT(1,2,3,line)
         do i = 1, NIND
            write(line,710) IVAR(M+i), XORIG(IVAR(M+i)), i
 710        format(' XORIG(',i5,') = ',d15.8,
     1             '  is independent number ',i5)
            call C_OUT(1,2,1,line)
         enddo
      endif
C
C     I think that's it...
C
 9999 continue
      return

 9008 IERR = 8
      call C_OUT(2,0,1,'Error while trying to open BASIS.DAT')
      goto 9999
 9009 IERR = 8
      call C_OUT(2,0,1,'Error while reading BASIS.DAT')
      goto 9999
      end

C ==============================================================================
C
C     Work space demand computation
C
C ==============================================================================

      subroutine PARTITION_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)

      implicit none
      include 'IPOPT.INC'
      integer N, M, NLB, NUB, NZA, LRW, LIW
      double precision DAT(*)
      integer IDAT(*)

      LRW = 0
      LIW = 0

      if( M.ne.0 .and. N-M.ne.0 ) then
         call CONSTR_WS(N, M, NLB, NUB, NZA, LRW, LIW, DAT, IDAT)
      endif

      LIW = max(LIW, N)
      return
      end
