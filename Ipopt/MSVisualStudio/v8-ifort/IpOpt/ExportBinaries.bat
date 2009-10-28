@echo off

mkdir ..\..\..\..\include
rmdir ..\..\..\..\include /s /q
mkdir ..\..\..\..\include\coin
copy ..\..\..\src\LinAlg\IpVector.hpp				..\..\..\..\include\coin /Y
copy ..\..\..\MSVisualStudio\v8-ifort\IpOpt\config_ipopt.h      ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpAlgTypes.hpp                     ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpCachedResults.hpp                    ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpDebug.hpp                            ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpException.hpp                        ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpIpoptApplication.hpp             ..\..\..\..\include\coin /Y
copy ..\..\..\src\Algorithm\IpIpoptCalculatedQuantities.hpp     ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpJournalist.hpp                       ..\..\..\..\include\coin /Y
copy ..\..\..\src\LinAlg\IpMatrix.hpp                           ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpNLP.hpp                          ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpObserver.hpp                         ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpoptConfig.h                          ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpOptionsList.hpp                      ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpReferenced.hpp                       ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpRegOptions.hpp                       ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpReturnCodes.h                    ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpReturnCodes.hpp                  ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpReturnCodes.inc                  ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpReturnCodes_inc.h                ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpSmartPtr.hpp                         ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpSolveStatistics.hpp              ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpStdCInterface.h                  ..\..\..\..\include\coin /Y
copy ..\..\..\src\LinAlg\IpSymMatrix.hpp                        ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpTaggedObject.hpp                     ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpTimedTask.hpp                        ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpTNLP.hpp                         ..\..\..\..\include\coin /Y
copy ..\..\..\src\Interfaces\IpTNLPReducer.hpp                  ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpTypes.hpp                            ..\..\..\..\include\coin /Y
copy ..\..\..\src\Common\IpUtils.hpp                            ..\..\..\..\include\coin /Y

mkdir ..\..\..\..\lib
rmdir ..\..\..\..\lib /s /q
mkdir ..\..\..\..\lib\win32
copy ..\release\Ipopt*.dll ..\..\..\..\lib\win32 /y
copy ..\release\Ipopt*.lib ..\..\..\..\lib\win32 /y

mkdir ..\..\..\..\lib\x64
copy ..\x64\release\Ipopt*.dll ..\..\..\..\lib\x64 /y
copy ..\x64\release\Ipopt*.lib ..\..\..\..\lib\x64 /y

copy README-LIB.TXT ..\..\..\..\lib\README.TXT /y
