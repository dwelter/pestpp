REM set variables used in script
set nslaves=9
set pestpp_file=this

REM create directories for PEST++ and the slaves to run in
xcopy /e /q /y template master\
FOR /L %%i IN (1,1,%nslaves%) DO (
  xcopy /e /q /y template slave%%i\
)

REM start YAMR master
start /D"%CD%\master" .\pest++ %pestpp_file% /H :4005

REM start YAMR slaves
FOR /L %%i IN (1,1,%nslaves%) DO (
 echo %%i
 start /D"%CD%\slave%%i" .\pest++ /H localhost:4005
)

