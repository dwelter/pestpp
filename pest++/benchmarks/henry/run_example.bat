REM set variables used in script
set nslaves=6
set pestpp_file=this

REM create directories for PEST++ and the slaves to run in
xcopy /e /q /y template master\
REM FOR /L %%i IN (1,1,%nslaves%) DO (
REM  xcopy /e /q /y template slave%%i\
REM )

REM start YAMR master
start /D"%CD%\master" .\pest++ %pestpp_file% /H :4006

REM start YAMR slaves
REM FOR /L %%i IN (1,1,%nslaves%) DO (
REM echo %%i
REM  start /D"%CD%\slave%%i" .\pest++ /H localhost:4006
REM )

