REM set variables used in script
set nslaves=6
set pestpp_file=this

REM create directories for PEST++ and the slaves to run in
xcopy /e /q /y template master\
 FOR /L %%i IN (1,1,%nslaves%) DO (
  xcopy /e /q /y template slave%%i\
 )

REM start YAMR master
start /D"%CD%\master" beopest64 %pestpp_file% /H :4008

REM start YAMR slaves
 FOR /L %%i IN (1,1,%nslaves%) DO (
 echo %%i
  start /D"%CD%\slave%%i" beopest64 %pestpp_file% /H localhost:4008
 )

