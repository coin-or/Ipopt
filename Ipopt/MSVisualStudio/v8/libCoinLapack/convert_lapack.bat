@echo off
rem Copyright (C) 2007 International Business Machines and others.
rem All Rights Reserved.
rem This file is distributed under the Eclipse Public License.
rem
rem Author:   Andreas Waechter              IBM     2007-06-14

@rem $Id: convert_lapack.bat 117 2007-06-15 16:54:36Z andreasw $

f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dgetf2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dgetrf.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dgetrs.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlae2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlaev2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlanst.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlansy.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlapy2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlarf.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlarfb.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlarfg.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlarft.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlartg.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlascl.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlaset.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlasr.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlasrt.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlassq.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlaswp.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dlatrd.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dorg2l.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dorg2r.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dorgql.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dorgqr.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dorgtr.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dpotf2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dpotrf.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dpotrs.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dsteqr.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dsterf.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dsyev.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dsytd2.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\dsytrd.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\ieeeck.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\ilaenv.f
f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\iparmq.f

move ..\..\..\..\ThirdParty\Lapack\LAPACK\SRC\*.c .

f2c ..\..\..\..\ThirdParty\Lapack\LAPACK\INSTALL\dlamch.f

move ..\..\..\..\ThirdParty\Lapack\LAPACK\INSTALL\dlamch.c .

