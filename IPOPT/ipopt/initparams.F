C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
C*******************************************************************************
C
      subroutine INITPARAMS(IERR, NARGS, ARGS, CARGS)
C
C*******************************************************************************
C
C    $Id: initparams.F 674 2005-06-30 20:54:16Z andreasw $
C
C-------------------------------------------------------------------------------
C                                 Title
C-------------------------------------------------------------------------------
C
CT    Determine algorithmic parameters - either from given ARGS or from
CT    file PARAMS.DAT
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
CP   IERR      O    INT    Error flag: <>0 Error
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

#include <config_f.h>
C
C-------------------------------------------------------------------------------
C                             Parameter list
C-------------------------------------------------------------------------------
C
      integer IERR
      integer NARGS
      double precision ARGS(*)
      character*(*) CARGS(*)
C
C-------------------------------------------------------------------------------
C                            Local variables
C-------------------------------------------------------------------------------
C
      integer NPARAMS           ! maximal number of parameters that can be
                                ! specified
#ifndef INCLUDE_CC
      parameter( NPARAMS = 80 )
#else
      parameter( NPARAMS = 83 )
#endif

      character*20 CPARAMS(NPARAMS)

      double precision PVALUE(NPARAMS)
      logical PSET(NPARAMS)     ! set to true, if corresponding parameters in
                                ! CPARAMS has been assigned a value

      integer FIND_PARAM, READ_ENTRY

      integer i, j, cnr, in, retval, ifile
      logical ex, new, is_qp
      character*20 name, st
      character*80 line
      double precision dp

      integer DIDMSG
      save DIDMSG
      data DIDMSG/0/
      data CPARAMS /            ! names of parameters
     1     'IFILE',             ! if set to 1, write output file
     1     'ICNR',              ! Fortran channal number for write output
     1     'IPRINT',            ! determines extend of output into file
     1     'IOUTPUT',           ! determines width of iteration output line
     1     'IMAXITER',          ! maximal number of iterations
     1     'IMAXCPUSEC',        ! maximal number of CPU seconds
     1     'DTOL',              ! expected error tolerance
     1     'ISCALERR',          ! determines whether error estimate is scaled
     1     'IERROR',            ! determines how error tolerance is to be
     1                          ! computed
     1     'DSCALMULTMAX',      ! s_max in error computation formula
     1     'IMUINIT',           ! determines how MU etc are to be initialized
     1     'ISCALE',            ! nonzero if automatic scaling is to be done
     1     'DSCALECUT',         ! cut-off value for automatic scaling
     1     'DFSCALE',           ! scaling factor for objective function
     1     'IFULL',             ! full or reduced-space option
     1     'IKKTSCALE',         ! how scaling of full KKT matrix is to be done
     1     'ISYMSOLV',          ! solver for indefinite linear system
     1     'ICG',               ! whether CG, and if so how to precondition
     1     'DCGTOL',            ! if CG, this is (relative preconditioned)
                                ! tolerance
     1     'IMAXCGITER',        ! maximum number of CG iteration (if not NIND)
     1     'IQUASI',            ! what kind of quasi-Newton method - or exact
     1     'IINITB',            ! initialization of reduced Hessian
     1     'ICORRECT',          ! how to correct reduced Hessian if indefinite
     1     'IDAMP',             ! damping of cross term or not
     1     'DSR1TOL',           ! when to switch to from BFGS to SR1
     1     'DSKIPFACT',         ! for BFGS skipping rule
     1     'IMERIT',            ! type of line search
     1     'ICNRM',             ! type of norm to measure constraint violation
     1     'IALPHA',            ! how to treat primal and dual fraction to the
     1                          ! boundary rule
     1     'IRESTO',            ! type of filter restoration phase
     1     'DRESTOKKTRED',      ! factor for sufficient reduction in KKT system
     1                          ! based reestoration phase
     1     'DS_F',              ! exponent in filter switching rule (obj. func.)
     1     'DS_THETA',          ! exponent in filter switching rule (constrts.)
     1     'DGAMMA_F',          ! factor in filter margin (obj. func.)
     1     'DGAMMA_THETA',      ! factor in filter margin(constrts.)
     1     'DDELTA',            ! weighting factor in switching condition
     1     'DLS_SAFE',          ! when to enforce only increasing penalty
     1                          ! parameters
     1     'DNUMIN',            ! minimal (starting) value of penalty parameter
     1     'DRHO',              ! parameter in update rule for penalty
     1                          ! parameter
     1     'DWATCHTOL',         ! when to activate watch dog technique
     1     'ISOC',              ! second order correction
     1     'IITERBLOCKMAX',     ! when to delete filter
     1     'ISELBAS',           ! how to select basis partition
     1     'DMAXCOND',          ! when to repartition
     1     'ILAMBDA',           ! how to compute multiplier estimates
     1     'DMU0',              ! initial value of barrier parameter
     1     'DMULIN',            ! factor for linear decrease of MU
     1     'DMUSUPER',          ! exponent for superlinear decrease of MU
     1     'DERRSUPER',         ! exponent for superlinear decrease of ERR
     1     'DMUERRFAC',         ! factor between MU and barrier error tolerance
     1     'DMAXERR',           ! upper bound for barrier error tolerance
     1     'DTAU',              ! factor in fraction to the boundary rule
     1     'DMOVEBOUNDS',       ! relative perturbation of bounds
     1     'DBNDFRAC',          ! minimal distance of starting point from
     1                          ! boundary as fraction of bound difference
     1     'DBNDPUSH',          ! absolute minimal distance of starting point
     1                          ! from bounds
     1     'IHESSVECT',         ! how to compute Hessian-vector products
     1     'DVCORRECTFACT',     ! factor for maximal deviation of V from S^-1 MU
     1     'DPIVTOL',           ! initial pivot tolerance in linear equation
                                ! solver
     1     'DPIVTOLMAX',        ! maximal pivot tolerance in linear equation
                                ! solver
     1     'DTRONCGTOL',        ! tolerance for CG within TRON (restoration
     1                          ! phase)
     1     'DTRONGTOL',         ! tolerance for 'gtol' in TRON (restoration
     1                          ! phase)
     1     'ITRONHESS',         ! how to compute Hessian-vector products in TRON
     1     'ITRON2DERIVS',      ! whether Gauss-Newton approximation in TRON
     1     'DDEPCONDIAG',       ! value for diagonal element in lower right
                                ! corner of KKT matrix when constraints are
                                ! linearly dependent (factor)
     1     'DMUDEPEXP',         ! exponent on MU for lower right corner of
                                ! of KKT matrix when constraints are
                                ! linearly dependent
                                ! corner of KKT matrix when constraints are
                                ! linearly dependent
     1     'IBWARMSTART',       ! If set to 1, the initial Quasi-Newton matrix
                                ! is read from file and at the end written to
                                ! file
     1     'ILMLEN',            ! number of (s,y) stored for limited memory
                                ! quasi-Newton update
     1     'DFILLINFACT',       ! factor for computation of workspace demand
                                ! for linear solvers (related to fill-in)
     1     'DRESTOREDFACT',     ! factor of minimal reduction of constraint
                                ! violation (2-norm) before exiting restoration
                                ! phase
     1     'DPRECFACT',         ! factor multiplied with machine precision to
                                ! compute relaxation at round-off sensitive
                                ! comparisons (e.g. line search)
     1     'DALPHAMINFRAC',     ! safety factor for ALPHA_MIN in filter line
                                ! search
     1     'DCONSTRRELAX',      ! minimum value for the constraint relaxation
                                ! (SE corner in KKT matrix, full-space only)
     1     'DLAMINITMAX',       ! if maxnorm of least square
                                ! estimate larger than this, set LAMBDA to zero
     1     'DTHETAMAXFACT',     ! factor for theta_max
                                !   (upper bound on feasibility in filter)
     1     'DTHETAMINFACT',     ! factor for theta_min
                                !   (part of switching rule)
     1     'DINFMAXTOL',        ! maximal allowed (unscaled) dual infeasibility
                                ! for successful termination
     1     'DLINDAMP',          ! value of linear damping parameter in
                                ! objective function
     1     'IREFINEITER',       ! minimum number of iterative refinement steps
     1     'IS_QP',             ! If 1, default values for QPs and LPs are
                                ! chosen
#ifndef INCLUDE_CC
     1     'DCMAXTOL'/          ! maximal allowed (unscaled) constraint
                                ! violation for successful termination
#else
     1     'DCMAXTOL',          !
     1     'IMPEC_TRIGGER',     ! option of triggering step modification for 
                                ! MPECs
     1     'DMPEC_ETAFACT',     ! factor multiplying the modification
     1     'DMPEC_THRESH'/      ! threshold mu below which modification is 
                                ! triggered
#endif

C
C*******************************************************************************
C
C                           Executable Statements
C
C*******************************************************************************
C
      do i = 1, NPARAMS
         PSET(i) = .false.
      enddo

      if( NARGS.ge.0 ) then
         inquire(file = 'PARAMS.DAT', exist = ex)
      else
         ex = .false.
      endif

      do j = 1, abs(NARGS)
         i = FIND_PARAM(NPARAMS, CPARAMS, CARGS(j))
         if( i.eq.0 ) then
            write(line,110) CARGS(j), ARGS(j)
            call C_OUT(0,0,1,line)
         elseif( PSET(i) ) then
            write(line,120) CPARAMS(i)
            call C_OUT(0,0,1,line)
         else
            PSET(i)   = .true.
            PVALUE(i) = ARGS(j)
         endif
      enddo

      if( ex ) then
C     Read the PARAMS.DAT file
         cnr = 1
         new = .true.
         open(cnr, file = 'PARAMS.DAT', status = 'old', err = 999 )
C     Loop to read all entries
 100     retval = READ_ENTRY(cnr, new, 's', dp, in, name)
         if( retval.eq.1 ) goto 200  ! end of file reached
         if( retval.ne.0 ) goto 9008 ! problem reading next entry
         retval = READ_ENTRY(cnr, new, 'd', dp, in, st)
         if( retval.ne.0 ) goto 9008
C     Find entry in CPARAMS list
         i = FIND_PARAM(NPARAMS, CPARAMS, name)
         if( i.eq.0 ) then
            write(line,110) name, dp
 110        format('Parameter ',a20,' with value ',d23.17,
     1           ' unknown.')
            call C_OUT(0,0,1,line)
            goto 100
         endif
         if( PSET(i) ) then
            write(line,120) CPARAMS(i), dp
 120        format('Parameter ',a20,
     1           ' has multiple assignments - using ',d12.5, ' .')
            call C_OUT(0,0,1,line)
         endif
C     Everything Ok, put entry into table
         PSET(i)   = .true.
         PVALUE(i) = dp
         goto 100
 200     continue
         close(cnr)
      endif
C
C     Now check every parameter if it has been assigned - other use default
C

C     IFILE
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IFILE')
      if( .not.PSET(i) ) then
         ifile = 0
      else
         ifile = int(PVALUE(i))
      endif
      if( ifile.lt.0 .or. ifile.gt.1 ) goto 9004

C     QCNR
      if( ifile.eq.0 ) then
         QCNR = 0
      else
         i = FIND_PARAM(NPARAMS, CPARAMS, 'ICNR')
         if( .not.PSET(i) .or. ifile.eq.0 ) then
            QCNR = 42
         else
            QCNR = int(PVALUE(i))
         endif
      endif
      if( QCNR.lt.0 ) goto 9004

C     Open the IPOPT.OUT output file
      if( QCNR.gt.0 ) then
         call C_OUT(-1,0,0,line)
      endif

C     QPRINT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IPRINT')
      if( .not.PSET(i) ) then
         QPRINT = 1
      else
         QPRINT = int(PVALUE(i))
      endif
      if( QPRINT.lt.-1 ) goto 9004

      if( QCNR.gt.1 .and. QPRINT.ge.1 ) then
         call C_OUT(1,0,1,'============= Algorithmic parameters ========
     1======')
         call C_OUT(1,0,1,' ')
         write(line,8000) 'ICNR',QCNR
         call C_OUT(1,0,1,line)
         write(line,8000) 'IPRINT',QPRINT
         call C_OUT(1,0,1,line)
      endif

C     QOUTPUT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IOUTPUT')
      if( .not.PSET(i) ) then
         QOUTPUT = 0
      else
         QOUTPUT = int(PVALUE(i))
      endif
      write(line,8000) 'IOUTPUT',QOUTPUT
      call C_OUT(1,1,1,line)
      if( QOUTPUT.lt.0 .or. QOUTPUT.gt.1 ) goto 9004

C     QMAXITER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMAXITER')
      if( .not.PSET(i) ) then
         QMAXITER = 10000
      else
         QMAXITER = int(PVALUE(i))
      endif
      write(line,8000) 'IMAXITER',QMAXITER
      call C_OUT(1,1,1,line)
      if( QMAXITER.lt.0 ) goto 9004

C     QMAXCPUSEC
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMAXCPUSEC')
      if( .not.PSET(i) ) then
         QMAXCPUSEC = 999999999
      else
         QMAXCPUSEC = int(PVALUE(i))
      endif
      write(line,8000) 'IMAXCPUSEC',QMAXCPUSEC
      call C_OUT(1,1,1,line)
      if( QMAXCPUSEC.le.0 ) goto 9004

C     QTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTOL')
      if( .not.PSET(i) ) then
         QTOL = 1.d-8
      else
         QTOL = PVALUE(i)
      endif
      write(line,8010) 'DTOL',QTOL
      call C_OUT(1,1,1,line)
      if( QTOL.le.0.d0 ) goto 9004

C     is_qp
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IS_QP')
      if( PSET(i) .and. PVALUE(i).ne.0.d0 ) then
         is_qp = .true.
         write(line,8000) 'IS_QP',int(PVALUE(i))
      else
         is_qp = .false.
         write(line,8000) 'IS_QP',0
      endif

C     QINIT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMUINIT')
      if( .not.PSET(i) ) then
         QINIT = 1
      else
         QINIT = int(PVALUE(i))
      endif
      write(line,8000) 'IMUINIT',QINIT
      call C_OUT(1,1,1,line)
      if( QINIT.lt.-1 .or. QINIT.gt.5 ) goto 9004

C     QSCALE
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ISCALE')
      if( .not.PSET(i) ) then
         QSCALE = 2
      else
         QSCALE = int(PVALUE(i))
      endif
      write(line,8000) 'ISCALE',QSCALE
      call C_OUT(1,1,1,line)
      if( QSCALE.lt.-1 .or. QSCALE.gt.5 ) goto 9004

C     QSCALECUT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DSCALECUT')
      if( .not.PSET(i) ) then
         if( QSCALE.ge.3 ) then
            QSCALECUT = 1.d-8
         else
            QSCALECUT = 100.d0
         endif
      else
         QSCALECUT = PVALUE(i)
      endif
      write(line,8010) 'DSCALECUT', QSCALECUT
      call C_OUT(1,1,1,line)
      if( QSCALECUT.le.0.d0 ) goto 9004

C     QFSCALE
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DFSCALE')
      if( .not.PSET(i) ) then
         QFSCALE = 1.d0
      else
         QFSCALE = PVALUE(i)
      endif
      write(line,8010) 'DFSCALE', QFSCALE
      call C_OUT(1,1,1,line)

C     QFULL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IFULL')
      if( .not.PSET(i) ) then
         QFULL = 1
      else
         QFULL = int(PVALUE(i))
      endif
      write(line,8000) 'IFULL', QFULL
      call C_OUT(1,1,1,line)
      if( QFULL.ne.0 .and. QFULL.ne.1 ) goto 9004

C     QSYMSOLV
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ISYMSOLV')
      if( .not.PSET(i) ) then
#ifndef HAVE_MA27
#ifndef HAVE_MA57
         QSYMSOLV = 0
#else
         QSYMSOLV = 2
#endif
#else
         QSYMSOLV = 1
#endif
      else
         QSYMSOLV = int(PVALUE(i))
      endif
      write(line,8000) 'ISYMSOLV', QSYMSOLV
      call C_OUT(1,1,1,line)
      if( QSYMSOLV.lt.0 .or. QSYMSOLV.gt.6 )
     1     goto 9004

C     QKKTSCALE
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IKKTSCALE')
      if( .not.PSET(i) ) then
         QKKTSCALE = 1
#ifndef HAVE_MC30
#ifndef HAVE_MC29
#ifndef HAVE_MC19
         QKKTSCALE = 0
#endif
#endif
#endif
      else
         QKKTSCALE = int(PVALUE(i))
      endif
      write(line,8000) 'IKKTSCALE', QKKTSCALE
      call C_OUT(1,1,1,line)
      if( QKKTSCALE.lt.0 .or. QKKTSCALE.gt.2 ) goto 9004

C     QSCALERR
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ISCALERR')
      if( .not.PSET(i) ) then
         QSCALERR = 1
      else
         QSCALERR = int(PVALUE(i))
      endif
      write(line,8000) 'ISCALERR', QSCALERR
      call C_OUT(1,1,1,line)
      if( QSCALERR.lt.0 .or. QSCALERR.gt.1 ) goto 9004

C     QERROR
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IERROR')
      if( .not.PSET(i) ) then
         if( QFULL.eq.0 ) then
            if( QSCALERR.eq.1 ) then
               QERROR = 1
            else
               QERROR = -1
            endif
         else
            if( QSCALERR.eq.1 ) then
               QERROR = 2
            else
               QERROR = -2
            endif
         endif
      else
         QERROR = int(PVALUE(i))
      endif
      write(line,8000) 'IERROR', QERROR
      call C_OUT(1,1,1,line)
      if( QERROR.lt.-2 .or. QERROR.gt.2 .or.
     1     (QFULL.ne.1 .and. abs(QERROR).eq.2) ) goto 9004

C     QSCALMULTMAX
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DSCALMULTMAX')
      if( .not.PSET(i) ) then
         QSCALMULTMAX = 1.d2
      else
         QSCALMULTMAX = PVALUE(i)
      endif
      write(line,8010) 'DSCALMULTMAX',QSCALMULTMAX
      call C_OUT(1,1,1,line)
      if( QSCALMULTMAX.le.0.d0 ) goto 9004

C     QCG
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ICG')
      if( .not.PSET(i) .or. QFULL.eq.1 ) then
         QCG = 0
      else
         QCG = int(PVALUE(i))
      endif
      write(line,8000) 'ICG', QCG
      call C_OUT(1,1,1,line)
      if( QCG.lt.-3 .or. QCG.gt.4 ) goto 9004

C     QCGTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DCGTOL')
      if( .not.PSET(i) ) then
         QCGTOL = 1.d-2
      else
         QCGTOL = PVALUE(i)
      endif
      write(line,8010) 'DCGTOL',QCGTOL
      call C_OUT(1,1,1,line)
      if( QCGTOL.eq.0.d0 .or. abs(QCGTOL).ge.1.d0 ) goto 9004

C     QMAXCGITER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMAXCGITER')
      if( .not.PSET(i) ) then
         QMAXCGITER = 0
      else
         QMAXCGITER = int(PVALUE(i))
      endif
      write(line,8000) 'IMAXCGITER', QMAXCGITER
      call C_OUT(1,1,1,line)
      if( QMAXCGITER.lt.0 ) goto 9004

C     QQUASI
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IQUASI')
      if( .not.PSET(i) ) then
         QQUASI = 0
      else
         QQUASI = int(PVALUE(i))
      endif
      write(line,8000) 'IQUASI', QQUASI
      call C_OUT(1,1,1,line)
      if( QQUASI.lt.-6 .or. QQUASI.gt.6 .or.
     1     (abs(QQUASI).eq.6.and.QFULL.eq.0) ) goto 9004

C     QLMLEN
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ILMLEN')
      if( .not.PSET(i) ) then
         QLMLEN = 6
      else
         QLMLEN = int(PVALUE(i))
      endif
      write(line,8000) 'ILMLEN', QLMLEN
      call C_OUT(1,1,1,line)
      if( QQUASI.eq.6 .and. QLMLEN.le.0 ) goto 9004

C     QINITB
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IINITB')
      if( .not.PSET(i) ) then
         QINITB = 0
      else
         QINITB = int(PVALUE(i))
      endif
      write(line,8000) 'IINITB', QINITB
      call C_OUT(1,1,1,line)
      if( QINITB.lt.0 .or. QINITB.gt.5 ) goto 9004

C     QBWARMSTART
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IBWARMSTART')
      if( .not.PSET(i) ) then
         QBWARMSTART = 0
      else
         QBWARMSTART = int(PVALUE(i))
      endif
      write(line,8000) 'IBWARMSTART', QBWARMSTART
      call C_OUT(1,1,1,line)
      if( QBWARMSTART.gt.1 .or. QBWARMSTART.lt.0 ) goto 9004

C     QCORRECT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ICORRECT')
      if( .not.PSET(i) ) then
         QCORRECT = 3
      else
         QCORRECT = int(PVALUE(i))
      endif
      write(line,8000) 'ICORRECT', QCORRECT
      call C_OUT(1,1,1,line)
      if( abs(QCORRECT).gt.4 .or. QCORRECT.eq.0 ) goto 9004

C     QDAMP
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IDAMP')
      if( .not.PSET(i) ) then
         QDAMP = 0
      else
         QDAMP = int(PVALUE(i))
      endif
      write(line,8000) 'IDAMP', QDAMP
      call C_OUT(1,1,1,line)
      if( QDAMP.lt.0 .or. QDAMP.gt.1 ) goto 9004

C     QSR1TOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DSR1TOL')
      if( .not.PSET(i) ) then
         QSR1TOL = 1.d0
      else
         QSR1TOL = PVALUE(i)
      endif
      write(line,8010) 'DSR1TOL', QSR1TOL
      call C_OUT(1,1,1,line)
      if( QSR1TOL.lt.0.d0 ) goto 9004

C     QSKIP
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DSKIPFACT')
      if( .not.PSET(i) ) then
         QSKIP = 1000.d0
      else
         QSKIP = PVALUE(i)
      endif
      write(line,8010) 'DSKIPFACT', QSKIP
      call C_OUT(1,1,1,line)
      if( QSKIP.lt.0.d0 ) goto 9004

C     QMERIT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMERIT')
      if( .not.PSET(i) ) then
         QMERIT = 4
      else
         QMERIT = int(PVALUE(i))
      endif
      write(line,8000) 'IMERIT', QMERIT
      call C_OUT(1,1,1,line)
      if( QMERIT.lt.-2 .or. QMERIT.gt.5 ) goto 9004

C     QCNRM
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ICNRM')
      if( .not.PSET(i) ) then
         QCNRM = 1
      else
         QCNRM = int(PVALUE(i))
      endif
      write(line,8000) 'ICNRM', QCNRM
      call C_OUT(1,1,1,line)
      if( QCNRM.lt.1 .or. QCNRM.gt.3 ) goto 9004

C     QALPHA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IALPHA')
      if( .not.PSET(i) ) then
         QALPHA = 2
      else
         QALPHA = int(PVALUE(i))
      endif
      write(line,8000) 'IALPHA', QALPHA
      call C_OUT(1,1,1,line)
      if( QALPHA.lt.0 .or. QALPHA.gt.2 ) goto 9004

C     QRESTO
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IRESTO')
      if( .not.PSET(i) ) then
         if( QFULL.eq.1 ) then
            QRESTO = 1
         else
            QRESTO = 2
         endif
      else
         QRESTO = int(PVALUE(i))
      endif
      write(line,8000) 'IRESTO', QRESTO
      call C_OUT(1,1,1,line)
      if( QRESTO.ne.1 .and. QRESTO.ne.2 ) goto 9004

C     QRESTOKKTRED
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DRESTOKKTRED')
      if( .not.PSET(i) ) then
         QRESTOKKTRED = 0.999d0
      else
         QRESTOKKTRED = PVALUE(i)
      endif
      write(line,8010) 'DRESTOKKTRED', QRESTOKKTRED
      call C_OUT(1,1,1,line)
      if( QRESTOKKTRED.lt.0.d0 .or. QRESTOKKTRED.ge.1.d0 ) goto 9004

C     QS_F
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DS_F')
      if( .not.PSET(i) ) then
C         QS_F = 1.1d0
         QS_F = 2.3d0
      else
         QS_F = PVALUE(i)
      endif
      write(line,8010) 'DS_F', QS_F
      call C_OUT(1,1,1,line)
      if( QS_F.lt.1.d0 ) goto 9004

C     QS_THETA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DS_THETA')
      if( .not.PSET(i) ) then
C         QS_THETA = 2.3d0
         QS_THETA = 1.1d0
      else
         QS_THETA = PVALUE(i)
      endif
      write(line,8010) 'DS_THETA', QS_THETA
      call C_OUT(1,1,1,line)
C      if( QS_F.le.2*QS_THETA ) goto 9004
      if( QS_THETA.le.1.d0 ) goto 9004

C     QGAMMA_F
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DGAMMA_F')
      if( .not.PSET(i) ) then
         QGAMMA_F = 1.d-5
      else
         QGAMMA_F = PVALUE(i)
      endif
      write(line,8010) 'DGAMMA_F', QGAMMA_F
      call C_OUT(1,1,1,line)
      if( QGAMMA_F.le.0.d0 .or. QGAMMA_F.ge.1.d0 ) goto 9004

C     QGAMMA_THETA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DGAMMA_THETA')
      if( .not.PSET(i) ) then
         QGAMMA_THETA = 1.d-5
      else
         QGAMMA_THETA = PVALUE(i)
      endif
      write(line,8010) 'DGAMMA_THETA', QGAMMA_THETA
      call C_OUT(1,1,1,line)
      if( QGAMMA_THETA.le.0.d0 .or. QGAMMA_THETA.ge.1.d0 ) goto 9004

C     QDELTA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DDELTA')
      if( .not.PSET(i) ) then
         QDELTA = 1.d0
      else
         QDELTA = PVALUE(i)
      endif
      write(line,8010) 'DDELTA', QDELTA
      call C_OUT(1,1,1,line)
      if( QDELTA.le.0.d0 ) goto 9004

C     QLS_SAFE
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DLS_SAFE')
      if( .not.PSET(i) ) then
         QLS_SAFE = 1.d20
      else
         QLS_SAFE = PVALUE(i)
      endif
      write(line,8010) 'DLS_SAFE', QLS_SAFE
      call C_OUT(1,1,1,line)
      if( QLS_SAFE.lt.0.d0 ) goto 9004

C     QNUMIN
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DNUMIN')
      if( .not.PSET(i) ) then
         QNUMIN = 1.d-6
      else
         QNUMIN = PVALUE(i)
      endif
      write(line,8010) 'DNUMIN', QNUMIN
      call C_OUT(1,1,1,line)
      if( QNUMIN.le.0.d0 ) goto 9004

C     QRHO
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DRHO')
      if( .not.PSET(i) ) then
         QRHO = 1.d-1
      else
         QRHO = PVALUE(i)
      endif
      write(line,8010) 'DRHO', QRHO
      call C_OUT(1,1,1,line)
      if( QRHO.le.0.d0 .or. QRHO.ge.1.d0 ) goto 9004

C     QWATCHTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DWATCHTOL')
      if( .not.PSET(i) ) then
         QWATCHTOL = 0.d0
      else
         QWATCHTOL = PVALUE(i)
      endif
      write(line,8010) 'DWATCHTOL', QWATCHTOL
      call C_OUT(1,1,1,line)

C     QSOC
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ISOC')
      if( .not.PSET(i) ) then
         if( abs(QMERIT).ne.5 ) then
            QSOC = 4
         else
            QSOC = 0
         endif
      else
         QSOC = int(PVALUE(i))
      endif
      write(line,8000) 'ISOC',QSOC
      call C_OUT(1,1,1,line)
      if( QSOC.lt.0 ) goto 9004

C     QITERBLOCKMAX
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IITERBLOCKMAX')
      if( .not.PSET(i) ) then
         QITERBLOCKMAX = 4
      else
         QITERBLOCKMAX = int(PVALUE(i))
      endif
      write(line,8000) 'IITERBLOCKMAX',QITERBLOCKMAX
      call C_OUT(1,1,1,line)
      if( QITERBLOCKMAX.lt.0 ) goto 9004

C     QSELBAS
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ISELBAS')
      if( .not.PSET(i) ) then
         QSELBAS = 0
      else
         QSELBAS = int(PVALUE(i))
      endif
      write(line,8000) 'ISELBAS',QSELBAS
      call C_OUT(1,1,1,line)
      if( QSELBAS.lt.0 .or. QSELBAS.gt.4 ) goto 9004

C     QMAXCOND
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMAXCOND')
      if( .not.PSET(i) ) then
         if( QSELBAS.eq.0 ) then
            QMAXCOND = 5.d1
         else
            QMAXCOND = 1.d20
         endif
      else
         QMAXCOND = PVALUE(i)
      endif
      write(line,8010) 'DMAXCOND',QMAXCOND
      call C_OUT(1,1,1,line)
      if( QMAXCOND.le.1.d0 ) goto 9004

C     QLAMBDA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ILAMBDA')
      if( .not.PSET(i) ) then
         if( QFULL.eq.0 ) then
            QLAMBDA = 1
         else
            QLAMBDA = -2
         endif
      else
         QLAMBDA = int(PVALUE(i))
      endif
      write(line,8000) 'ILAMBDA', QLAMBDA
      call C_OUT(1,1,1,line)
      if( QLAMBDA.lt.-2 .or. QLAMBDA.gt.2 .or. QLAMBDA.eq.-1 .or.
     1     (QFULL.ne.1 .and. abs(QLAMBDA).eq.2) ) goto 9004

C     QMUINIT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMU0')
      if( .not.PSET(i) ) then
         if( is_qp ) then
            QMUINIT = 1.d3
         else
            QMUINIT = 1.d-1
         endif
      else
         QMUINIT = PVALUE(i)
      endif
      write(line,8010) 'DMU0',QMUINIT
      call C_OUT(1,1,1,line)
      if( QMUINIT.le.0.d0 ) goto 9004

C     QMULIN
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMULIN')
      if( .not.PSET(i) ) then
         QMULIN = 2.d-1
      else
         QMULIN = PVALUE(i)
      endif
      write(line,8010) 'DMULIN',QMULIN
      call C_OUT(1,1,1,line)
      if( QMULIN.le.0.d0 .or. QMULIN.ge.1.d0 ) goto 9004

C     QMUSUPER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMUSUPER')
      if( .not.PSET(i) ) then
         QMUSUPER = 1.5d0
      else
         QMUSUPER = PVALUE(i)
      endif
      write(line,8010) 'DMUSUPER',QMUSUPER
      call C_OUT(1,1,1,line)
      if( QMUSUPER.lt.1.d0 ) goto 9004

C     QERRSUPER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DERRSUPER')
      if( .not.PSET(i) ) then
         QERRSUPER = 1.1d0
      else
         QERRSUPER = PVALUE(i)
      endif
      write(line,8010) 'DERRSUPER',QERRSUPER
      call C_OUT(1,1,1,line)
      if( QERRSUPER.lt.1.d0 ) goto 9004

C     QMUERRFAC
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMUERRFAC')
      if( .not.PSET(i) ) then
         QMUERRFAC = 1.d1
      else
         QMUERRFAC = PVALUE(i)
      endif
      write(line,8010) 'DMUERRFAC',QMUERRFAC
      call C_OUT(1,1,1,line)
      if( QMUERRFAC.le.0.d0 ) goto 9004

C     QMAXERR
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMAXERR')
      if( .not.PSET(i) ) then
         QMAXERR = 1.d300
      else
         QMAXERR = PVALUE(i)
      endif
      write(line,8010) 'DMAXERR',QMAXERR
      call C_OUT(1,1,1,line)
      if( QMAXERR.lt.QTOL ) goto 9004

C     QTAU
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTAU')
      if( .not.PSET(i) ) then
         QTAU = 0.99d0
      else
         QTAU = PVALUE(i)
      endif
      write(line,8010) 'DTAU',QTAU
      QTAUMAX = QTAU
      call C_OUT(1,1,1,line)
      if( QTAU.le.0.d0 .or. QTAU.ge.1.d0 ) goto 9004

C     QMOVEBOUNDS
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMOVEBOUNDS')
      if( .not.PSET(i) ) then
         QMOVEBOUNDS = 1.d-8
      else
         QMOVEBOUNDS = PVALUE(i)
      endif
      write(line,8010) 'DMOVEBOUNDS',QMOVEBOUNDS
      call C_OUT(1,1,1,line)
      if( QMOVEBOUNDS.lt.0.d0 ) goto 9004

C     QFRAC
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DBNDFRAC')
      if( .not.PSET(i) ) then
         if( is_qp ) then
            QFRAC = 1.d-1
         else
            QFRAC = 1.d-2
         endif
      else
         QFRAC = PVALUE(i)
      endif
      write(line,8010) 'DBNDFRAC',QFRAC
      call C_OUT(1,1,1,line)
      if( QFRAC.le.0.d0 .or. QFRAC.ge.5.d-1 ) goto 9004

C     QTHETA
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DBNDPUSH')
      if( .not.PSET(i) ) then
         if( is_qp ) then
            QTHETA = 1.d-1
         else
            QTHETA = 1.d-2
         endif
      else
         QTHETA = PVALUE(i)
      endif
      write(line,8010) 'DBNDPUSH',QTHETA
      call C_OUT(1,1,1,line)
      if( QTHETA.le.0.d0 ) goto 9004

C     QHESS
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IHESSVECT')
      if( .not.PSET(i) ) then
         if( QFULL.eq.0 ) then
            QHESS = 1
         else
            QHESS = 0
         endif
      else
         QHESS = int(PVALUE(i))
      endif
      write(line,8000) 'IHESSVECT',QHESS
      call C_OUT(1,1,1,line)
      if( QHESS.lt.0 .or. QHESS.gt.2 ) goto 9004

C     QVCORRECTFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DVCORRECTFACT')
      if( .not.PSET(i) ) then
         QVCORRECTFACT = 1.d10
      else
         QVCORRECTFACT = PVALUE(i)
      endif
      write(line,8010) 'DVCORRECTFACT',QVCORRECTFACT
      call C_OUT(1,1,1,line)
      if( QVCORRECTFACT.lt.1.d0 ) goto 9004

C     QPIVTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DPIVTOL')
      if( .not.PSET(i) ) then
         QPIVTOL = 1.d-8
      else
         QPIVTOL = PVALUE(i)
      endif
      write(line,8010) 'DPIVTOL',QPIVTOL
      call C_OUT(1,1,1,line)
      if( QPIVTOL.le.0.d0 .or. QPIVTOL.gt.5.d-1 ) goto 9004

C     QPIVTOLMAX
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DPIVTOLMAX')
      if( .not.PSET(i) ) then
         QPIVTOLMAX = 1.d-6
      else
         QPIVTOLMAX = PVALUE(i)
      endif
      if( QSYMSOLV.eq.3 .or. QSYMSOLV.eq.4 ) then
         QPIVTOLMAX = QPIVTOL   ! There is no PIVTOL for WSMP
      endif
      if( QPIVTOLMAX.lt.QPIVTOL ) then
         write(line,*)
     1 'Warning: Increasing QPIVTOLMAX to value of QPIVTOL = ',QPIVTOL
         call C_OUT(1,0,1,line)
         QPIVTOLMAX = QPIVTOL
      endif
      write(line,8010) 'DPIVTOLMAX',QPIVTOLMAX
      call C_OUT(1,1,1,line)
      if( QPIVTOLMAX.lt.QPIVTOL .or. QPIVTOLMAX.gt.5.d-1 ) goto 9004

C     QTRONCGTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTRONCGTOL')
      if( .not.PSET(i) ) then
         QTRONCGTOL = 1.d-3
      else
         QTRONCGTOL = PVALUE(i)
      endif
      write(line,8010) 'DTRONCGTOL',QTRONCGTOL
      call C_OUT(1,1,1,line)
      if( QTRONCGTOL.le.0.d0 .or. QTRONCGTOL.ge.1.d0 ) goto 9004

C     QTRONGTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTRONGTOL')
      if( .not.PSET(i) ) then
         QTRONGTOL = 1.d-3*QTOL
      else
         QTRONGTOL = PVALUE(i)
      endif
      write(line,8010) 'DTRONGTOL',QTRONGTOL
      call C_OUT(1,1,1,line)
      if( QTRONGTOL.le.0.d0 ) goto 9004

C     QTRONHESS
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ITRONHESS')
      if( .not.PSET(i) ) then
         if( QFULL.eq.0 ) then
            QTRONHESS = 0
         else
            QTRONHESS = 0
#ifdef HAVE_MC35
#ifdef HAVE_MC39
            QTRONHESS = 1
#endif
#endif
         endif
      else
         QTRONHESS = int(PVALUE(i))
      endif
      write(line,8000) 'ITRONHESS',QTRONHESS
      call C_OUT(1,1,1,line)
      if( QTRONHESS.lt.0 .or. QTRONHESS.gt.2 ) goto 9004

C     QTRON2DERIVS
      i = FIND_PARAM(NPARAMS, CPARAMS, 'ITRON2DERIVS')
      if( .not.PSET(i) ) then
         QTRON2DERIVS = 1
      else
         QTRON2DERIVS = int(PVALUE(i))
      endif
      write(line,8000) 'ITRON2DERIVS',QTRON2DERIVS
      call C_OUT(1,1,1,line)
      if( QTRON2DERIVS.lt.0 .or. QTRON2DERIVS.gt.1 ) goto 9004

C     QDEPCONDIAG
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DDEPCONDIAG')
      if( .not.PSET(i) ) then
         QDEPCONDIAG = 1.d-8
      else
         QDEPCONDIAG = PVALUE(i)
      endif
      write(line,8010) 'DDEPCONDIAG',QDEPCONDIAG
      call C_OUT(1,1,1,line)
      if( QDEPCONDIAG.lt.0.d0 ) goto 9004

C     QMUDEPEXP
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMUDEPEXP')
      if( .not.PSET(i) ) then
         QMUDEPEXP = 0.25d0
      else
         QMUDEPEXP = PVALUE(i)
      endif
      write(line,8010) 'DMUDEPEXP',QMUDEPEXP
      call C_OUT(1,1,1,line)
      if( QMUDEPEXP.lt.0.d0 ) goto 9004

C     QFILLINFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DFILLINFACT')
      if( .not.PSET(i) ) then
         if( QSYMSOLV.le.2 ) then
            QFILLINFACT = 5.d0
         else
            QFILLINFACT = 20.d0
         endif
      else
         QFILLINFACT = PVALUE(i)
      endif
      write(line,8010) 'DFILLINFACT',QFILLINFACT
      call C_OUT(1,1,1,line)
      if( QFILLINFACT.lt.1.d0 ) goto 9004

C     QRESTOREDFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DRESTOREDFACT')
      if( .not.PSET(i) ) then
         if( is_qp ) then
            QRESTOREDFACT = 1.d-1
         else
            QRESTOREDFACT = 9.d-1
         endif
      else
         QRESTOREDFACT = PVALUE(i)
      endif
      write(line,8010) 'DRESTOREDFACT',QRESTOREDFACT
      call C_OUT(1,1,1,line)
      if( QRESTOREDFACT.lt.0.d0 ) goto 9004

C     QPRECFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DPRECFACT')
      if( .not.PSET(i) ) then
         QPRECFACT = 10.d0
      else
         QPRECFACT = PVALUE(i)
      endif
      write(line,8010) 'DPRECFACT',QPRECFACT
      call C_OUT(1,1,1,line)
      if( QPRECFACT.lt.0.d0 ) goto 9004

C     QALPHAMINFRAC
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DALPHAMINFRAC')
      if( .not.PSET(i) ) then
         QALPHAMINFRAC = 5.d-2
      else
         QALPHAMINFRAC = PVALUE(i)
      endif
      write(line,8010) 'DALPHAMINFRAC',QALPHAMINFRAC
      call C_OUT(1,1,1,line)
      if( QALPHAMINFRAC.le.0.d0 ) goto 9004

C     QCONSTRRELAX
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DCONSTRRELAX')
      if( .not.PSET(i) ) then
         QCONSTRRELAX = 0.d0
      else
         QCONSTRRELAX = PVALUE(i)
      endif
      write(line,8010) 'DCONSTRRELAX',QCONSTRRELAX
      call C_OUT(1,1,1,line)
      if( QCONSTRRELAX.lt.0.d0 ) goto 9004

C     QLAMINITMAX
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DLAMINITMAX')
      if( .not.PSET(i) ) then
         QLAMINITMAX = 1.d3
      else
         QLAMINITMAX = PVALUE(i)
      endif
      write(line,8010) 'DLAMINITMAX',QLAMINITMAX
      call C_OUT(1,1,1,line)
      if( QLAMINITMAX.lt.0.d0 ) goto 9004

C     QLAMINIT
      if( QLAMINITMAX.eq.0.d0 ) then
         QLAMINIT = 0
      else
         QLAMINIT = 1
      endif

C     QREFINEITER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IREFINEITER')
      if( .not.PSET(i) ) then
         QREFINEITER = 1
      else
         QREFINEITER = int(PVALUE(i))
      endif
      write(line,8000) 'IREFINEITER',QREFINEITER
      call C_OUT(1,1,1,line)
      if( QREFINEITER.lt.0 .) goto 9004

C     QTHETAMAXFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTHETAMAXFACT')
      if( .not.PSET(i) ) then
         QTHETAMAXFACT = 1.d4
      else
         QTHETAMAXFACT = PVALUE(i)
      endif
      write(line,8010) 'DTHETAMAXFACT',QTHETAMAXFACT
      call C_OUT(1,1,1,line)
      if( QTHETAMAXFACT.lt.1.d0 ) goto 9004

C     QTHETAMINFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DTHETAMINFACT')
      if( .not.PSET(i) ) then
         QTHETAMINFACT = 1.d-4
      else
         QTHETAMINFACT = PVALUE(i)
      endif
      write(line,8010) 'DTHETAMINFACT',QTHETAMINFACT
      call C_OUT(1,1,1,line)
      if( QTHETAMINFACT.le.0.d0 ) goto 9004

C     QDINFMAXTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DINFMAXTOL')
      if( .not.PSET(i) ) then
C         QDINFMAXTOL = sqrt(QTOL)
         QDINFMAXTOL = 1.d300
      else
         QDINFMAXTOL = PVALUE(i)
      endif
      write(line,8010) 'DINFMAXTOL',QDINFMAXTOL
      call C_OUT(1,1,1,line)
      if( QDINFMAXTOL.le.0.d0 ) goto 9004

C     QCMAXTOL
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DCMAXTOL')
      if( .not.PSET(i) ) then
C         QCMAXTOL = sqrt(QTOL)
         QCMAXTOL = 1.d300
      else
         QCMAXTOL = PVALUE(i)
      endif
      write(line,8010) 'DCMAXTOL',QCMAXTOL
      call C_OUT(1,1,1,line)
      if( QCMAXTOL.le.0.d0 ) goto 9004

C     QLINDAMP
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DLINDAMP')
      if( .not.PSET(i) ) then
         if( QQUASI.eq.0 ) then
            QLINDAMP = 1.d-4
         else
            QLINDAMP = 0.d0
         endif
      else
         QLINDAMP = PVALUE(i)
      endif
      write(line,8010) 'DLINDAMP',QLINDAMP
      call C_OUT(1,1,1,line)
      if( QLINDAMP.lt.0.d0 ) goto 9004

#ifdef INCLUDE_CC
C     QMPEC_TRIGGER
      i = FIND_PARAM(NPARAMS, CPARAMS, 'IMPEC_TRIGGER')
      if( .not.PSET(i) ) then
         QMPEC_TRIGGER = 0
      else
         QMPEC_TRIGGER = PVALUE(i)
      endif
      write(line,8000) 'IMPEC_TRIGGER',QMPEC_TRIGGER
      call C_OUT(1,1,1,line)
      if( QMPEC_TRIGGER.lt.0 ) goto 9004

C     QMPEC_ETAFACT
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMPEC_ETAFACT')
      if( .not.PSET(i) ) then
         QMPEC_ETAFACT = 0.1d0
      else
         QMPEC_ETAFACT = PVALUE(i)
      endif
      write(line,8010) 'DMPEC_ETAFACT',QMPEC_ETAFACT
      call C_OUT(1,0,1,line)
      if( QMPEC_ETAFACT.lt.0.d0 ) goto 9004

C     QMPEC_THRESH
      i = FIND_PARAM(NPARAMS, CPARAMS, 'DMPEC_THRESH')
      if( .not.PSET(i) ) then
         QMPEC_THRESH = 5.d-6
      else
         QMPEC_THRESH = PVALUE(i)
      endif
      write(line,8010) 'DMPEC_THRESH',QMPEC_THRESH
      call C_OUT(1,0,1,line)
      if( QMPEC_THRESH.lt.0.d0 ) goto 9004
#else
      QMPEC_TRIGGER = 0
      QMPEC_ETAFACT = 0.d0
      QMPEC_THRESH  = 0.d0
#endif

      if( QCNR.gt.1 .and. QPRINT.ge.1 ) then
         call C_OUT(1,0,1,' ')
         call C_OUT(1,0,1,'========= End of algorithmic parameters =====
     1======')
         call C_OUT(1,0,1,' ')
      endif

C     Necessary for logic in algorithm
      if( QMERIT.eq.0 ) QSOC = 0

      if( QFULL.eq.1 .and. QMERIT.eq.2 .and. QQUASI.eq.0 ) then
         QKNITROLS = 1
         inquire(file='LINEAR',exist=ex)
         if( ex ) then
            QKNITROLIN = 1
            call C_OUT(2,0,1,'LINEAR MODEL IN LINE SEARCH!')
         else
            QKNITROLIN = 0
         endif
      else
         QKNITROLS = 0
      endif

      if( DIDMSG.eq.0 ) then
        DIDMSG = 1
        call C_OUT(2,-1000,1,'******************************************
     1************************************')
        call C_OUT(2,-1000,1,
     1'This program contains IPOPT, a program for large-scale nonlinear 
     2optimization.')
        call C_OUT(2,-1000,1,
     1'   IPOPT is released as open source under the Common Public Licen
     2se (CPL).')
        call C_OUT(2,-1000,1,
     1'               For more information visit www.coin-or.org/Ipopt')
        call C_OUT(2,-1000,1,'******************************************
     1************************************')
      endif

C     Print version number
      if( QPRINT.ge.1) then
         call C_OUT(1,1,1,'This is IPOPT version 2.2.1e')
      endif

#ifdef HAVE_VERSION_INC
C     Print out compilation information if desired
      if( QPRINT.ge.2 .and. QCNR.gt.0 ) then
         call C_OUT(1,0,1,' ')
         call C_OUT(1,0,1,'============= Compilation information =======
     1======')
         call C_OUT(1,0,1,' ')
         include 'VERSION.INC'
         call C_OUT(1,0,1,' ')
         call C_OUT(1,0,1,'========== End of compilation information ===
     1======')
         call C_OUT(1,0,1,' ')
      endif
#endif

 8000 format(2x,a20,' = ',i23)
 8010 format(2x,a20,' = ',d23.17)

      IERR = 0

C     Initialize a few common block variables
      FEVALS = 0
      CEVALS = 0
      LASTPRINT = -1

      return

 998  continue
      call C_OUT(0,0,1,'Error: Can''t open IPOPT.OUT.')
      IERR = 8
      return

 999  continue
      call C_OUT(0,0,1,'Error: Can''t open PARAMS.DAT.')
      IERR = 8
      return

 9004 IERR = 4
      write(line,*) 'initparams: Problems with IPOPT parameter ',
     1     CPARAMS(i)
      call C_OUT(2,0,1,line)
      return

 9008 IERR = 8
      return

      end
