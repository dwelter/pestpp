REM set variables used in script
set nslaves=10
set host=localhost
set port=4005
set pestpp_file=storage5

REM create directories for PEST++ and the gslave to run in
#xcopy /e /q /y template master\
FOR /L %%i IN (1,1,%nslaves%) DO (
  xcopy /e /q /y template slave%%i\
)

#REM start YAMR master
#start /D"%CD%\master" .\pest++ %pestpp_file% /H :%port%

REM start YAMR slaves
FOR /L %%i IN (1,1,%nslaves%) DO (
 echo %%i
 start /D"%CD%\slave%%i" .\pest++ /H %host%:%port%
)

