@echo off
rem Copyright (C) 2007 International Business Machines and others.
rem All Rights Reserved.
rem This file is distributed under the Eclipse Public License.
rem
rem Author:   Andreas Waechter              IBM     2007-06-14

@rem $Id: convert_hsl.bat 117 2007-06-15 16:54:36Z andreasw $

f2c ..\..\..\..\ThirdParty\HSL\ma27ad.f
f2c ..\..\..\..\ThirdParty\HSL\mc19ad.f

move ..\..\..\..\ThirdParty\HSL\*.c .
