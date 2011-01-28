@echo off
rem Copyright (C) 2007 International Business Machines and others.
rem All Rights Reserved.
rem This file is distributed under the Eclipse Public License.
rem
rem Author:   Andreas Waechter              IBM     2007-06-14

@rem $Id: convert_blas.bat 117 2007-06-15 16:54:36Z andreasw $

f2c ..\..\..\..\ThirdParty\Blas\dasum.f
f2c ..\..\..\..\ThirdParty\Blas\daxpy.f
f2c ..\..\..\..\ThirdParty\Blas\dcopy.f
f2c ..\..\..\..\ThirdParty\Blas\ddot.f
f2c ..\..\..\..\ThirdParty\Blas\dgbmv.f
f2c ..\..\..\..\ThirdParty\Blas\dgemm.f
f2c ..\..\..\..\ThirdParty\Blas\dgemv.f
f2c ..\..\..\..\ThirdParty\Blas\dger.f
f2c ..\..\..\..\ThirdParty\Blas\dnrm2.f
f2c ..\..\..\..\ThirdParty\Blas\drot.f
f2c ..\..\..\..\ThirdParty\Blas\drotg.f
f2c ..\..\..\..\ThirdParty\Blas\drotm.f
f2c ..\..\..\..\ThirdParty\Blas\drotmg.f
f2c ..\..\..\..\ThirdParty\Blas\dsbmv.f
f2c ..\..\..\..\ThirdParty\Blas\dscal.f
f2c ..\..\..\..\ThirdParty\Blas\dsdot.f
f2c ..\..\..\..\ThirdParty\Blas\dspmv.f
f2c ..\..\..\..\ThirdParty\Blas\dspr2.f
f2c ..\..\..\..\ThirdParty\Blas\dspr.f
f2c ..\..\..\..\ThirdParty\Blas\dswap.f
f2c ..\..\..\..\ThirdParty\Blas\dsymm.f
f2c ..\..\..\..\ThirdParty\Blas\dsymv.f
f2c ..\..\..\..\ThirdParty\Blas\dsyr2.f
f2c ..\..\..\..\ThirdParty\Blas\dsyr2k.f
f2c ..\..\..\..\ThirdParty\Blas\dsyr.f
f2c ..\..\..\..\ThirdParty\Blas\dsyrk.f
f2c ..\..\..\..\ThirdParty\Blas\dtbmv.f
f2c ..\..\..\..\ThirdParty\Blas\dtbsv.f
f2c ..\..\..\..\ThirdParty\Blas\dtpmv.f
f2c ..\..\..\..\ThirdParty\Blas\dtpsv.f
f2c ..\..\..\..\ThirdParty\Blas\dtrmm.f
f2c ..\..\..\..\ThirdParty\Blas\dtrmv.f
f2c ..\..\..\..\ThirdParty\Blas\dtrsm.f
f2c ..\..\..\..\ThirdParty\Blas\dtrsv.f
f2c ..\..\..\..\ThirdParty\Blas\idamax.f
f2c ..\..\..\..\ThirdParty\Blas\lsame.f
f2c ..\..\..\..\ThirdParty\Blas\xerbla.f

move ..\..\..\..\ThirdParty\Blas\*.c .
