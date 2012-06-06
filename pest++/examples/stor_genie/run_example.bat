REM set variables used in script
set nslaves=3
set gman_ip=127.0.0.1
set gman_port=24772
set pestpp_file=storage5

REM create directories for PEST++ and the gslave to run in
xcopy /q /y template master\
FOR /L %%i IN (1,1,%nslaves%) DO (
  xcopy /q /y template slave%%i\
)

REM start GMAN
start  /D"%CD%" gman /ip %gman_ip% /port %gman_port%

REM start GSLAVES
FOR /L %%i IN (1,1,%nslaves%) DO (
 echo %%i
  start /D"%CD%\slave%%i" gslave /interval 1.0 /console on /name Slave%%i /host %gman_ip%:%gman_port%
)

REM start PEST++
start /D"%CD%\master" pest++ %pestpp_file%
