C Copyright (C) 2002, Carnegie Mellon University and others.
C All Rights Reserved.
C This code is published under the Common Public License.
      double precision TIME_BB, TIME_CG, TIME_YPY, TIME_EXACTW
      double precision TIME_ZWZY_BACKS, TIME_ZWZY_EVALA
      double precision TIME_PZ_CHOL, TIME_HV, TIME_GET_STEP_FULL
      integer COUNT_CG, COUNT_RESTO_ITER, COUNT_NEG_CURV
      integer COUNT_RESTO_CALL, COUNT_TRON_CG, COUNT_HV
      integer COUNT_DEPCON
      common /CPUTIMES/ TIME_BB, TIME_CG, TIME_YPY, TIME_EXACTW,
     1     TIME_ZWZY_BACKS, TIME_ZWZY_EVALA, TIME_PZ_CHOL, TIME_HV,
     2     TIME_GET_STEP_FULL,
     2     COUNT_CG, COUNT_RESTO_ITER, COUNT_NEG_CURV,
     3     COUNT_RESTO_CALL, COUNT_TRON_CG, COUNT_HV, COUNT_DEPCON
      save   /CPUTIMES/
