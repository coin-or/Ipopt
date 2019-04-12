C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.

C $Id: IPOPT.INC 632 2004-07-20 23:52:15Z andreasw $

      integer          QMAXITER
      integer          QMAXCPUSEC
      integer          QCNR
      integer          QPRINT
      integer          QOUTPUT
      integer          QSELBAS
      integer          QQUASI
      integer          QINITB
      integer          QMERIT
      integer          QCNRM
      integer          QCORRECT
      integer          QFULL
      integer          QLAMBDA
      integer          QERROR
      integer          QSCALERR
      integer          QSOC
      integer          QRESTO
      integer          QINIT
      integer          QSCALE
      integer          QCG
      integer          QALPHA
      integer          QHESS
      integer          QDAMP
      integer          QTRONHESS
      integer          QTRON2DERIVS
      integer          QITERBLOCKMAX
      integer          QKKTSCALE
      integer          QSYMSOLV
      integer          QMAXCGITER
      integer          QBWARMSTART
      integer          QLMLEN
      integer          QLAMINIT
      integer          INMEMCHECK
      integer          QMPEC_TRIGGER
      integer          QREFINEITER

      double precision QTOL
      double precision QWATCHTOL
      double precision QLS_SAFE
      double precision QRHO
      double precision QNUMIN
      double precision QMAXCOND
      double precision QSR1TOL
      double precision QMUINIT
      double precision QMULIN
      double precision QMUSUPER
      double precision QMAXERR
      double precision QMUERRFAC
      double precision QMOVEBOUNDS
      double precision QFRAC
      double precision QTHETA
      double precision QTAU
      double precision QFSCALE
      double precision QSKIP
      double precision QVCORRECTFACT
      double precision QPIVTOL
      double precision QPIVTOLMAX
      double precision QTRONCGTOL
      double precision QTRONGTOL
      double precision QERRSUPER
      double precision QDEPCONDIAG
      double precision QRESTOKKTRED
      double precision QS_F
      double precision QGAMMA_F
      double precision QS_THETA
      double precision QGAMMA_THETA
      double precision QDELTA
      double precision QCGTOL
      double precision QFILLINFACT
      double precision QRESTOREDFACT
      double precision QPRECFACT
      double precision QALPHAMINFRAC
      double precision QCONSTRRELAX
      double precision QDINFMAXTOL
      double precision QCMAXTOL
      double precision QTHETAMAXFACT
      double precision QTHETAMINFACT
      double precision QLINDAMP
      double precision QMUDEPEXP
      double precision QLAMINITMAX
      double precision QSCALMULTMAX
      double precision QSCALECUT
      double precision QTAUMAX

      double precision QDDMAX
      parameter      ( QDDMAX = 1d0 )

      double precision QMPEC_ETAFACT
      double precision QMPEC_THRESH

      common /PARAMS/
     1     QTOL, QWATCHTOL, QLS_SAFE, QRHO, QNUMIN, QMAXCOND,
     1     QSR1TOL, QMUINIT, QMULIN, QMUSUPER, QMAXERR, QMUERRFAC,
     1     QMOVEBOUNDS,
     1     QFRAC, QTHETA, QTAU, QFSCALE, QSKIP, QVCORRECTFACT,
     1     QPIVTOL, QTRONCGTOL, QERRSUPER, QDEPCONDIAG, QRESTOKKTRED,
     1     QS_F, QGAMMA_F, QS_THETA, QGAMMA_THETA, QTRONGTOL, QCGTOL,
     1     QPIVTOLMAX, QFILLINFACT, QRESTOREDFACT, QPRECFACT,
     1     QALPHAMINFRAC, QCONSTRRELAX, QDINFMAXTOL, QCMAXTOL,
     1     QTHETAMAXFACT, QTHETAMINFACT, QLINDAMP, QMUDEPEXP,
     1     QLAMINITMAX, QSCALMULTMAX, QSCALECUT, QDELTA, QTAUMAX,
     1     QMAXITER, QCNR, QPRINT, QOUTPUT, QSELBAS, QQUASI, QINITB,
     1     QMERIT, QCNRM, QCORRECT, QFULL, QLAMBDA,
     1     QERROR, QSCALERR, QSOC, QRESTO, QINIT, QSCALE,
     1     QCG, QALPHA, QHESS, QDAMP, QTRONHESS, QTRON2DERIVS,
     1     QITERBLOCKMAX, QKKTSCALE, QSYMSOLV, QMAXCGITER, QBWARMSTART,
     1     QLMLEN, QLAMINIT, QREFINEITER, QMAXCPUSEC, INMEMCHECK
      save /PARAMS/
C
C     For AMPL complementarity constraints
C
      double precision AMPLMU
      common /AMPLMU/ AMPLMU
      save /AMPLMU/

C     MPEC parameters
      common /MPECPARAMS/ QMPEC_ETAFACT, QMPEC_THRESH, QMPEC_TRIGGER
      save /MPECPARAMS/

C
C     Constants for indefiniteness regularization
C
      double precision REGU_INIT, REGU_MAX, REGU_DEC_FACT, REGU_INC_FACT
      double precision REGU_INIT_FACT, REGU_MIN
      parameter( REGU_INIT = 1.d-4, REGU_MAX = 1.d40,
     1     REGU_DEC_FACT = 1.d0/3.d0, REGU_INC_FACT = 8.d0,
     1     REGU_INIT_FACT = 100.d0, REGU_MIN = 1.d-20 )

C
C     Set the following to true for memory debugging
C
      logical MEMDBG
      parameter( MEMDBG = .false. )
C
C     Set the following to true for some timing and memory information
C
      logical PRINTTIMING
      parameter( PRINTTIMING = .false. )

C
C     Following is set to one if we are performing a KNITRO-DIRECT type
C     of linesearch
C
      integer QKNITROLS, QKNITROLIN
      common /knitrols/ QKNITROLS, QKNITROLIN
      save /knitrols/

      double precision QLINDAMPEXP
      parameter( QLINDAMPEXP = 1.d0 )
c      parameter( QLINDAMPEXP = 5.d-1 )

C     Counting number of function evaluations
      integer FEVALS, CEVALS, LASTPRINT
      common /EVALS/ FEVALS, CEVALS, LASTPRINT
      save /EVALS/
