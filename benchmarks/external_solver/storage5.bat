set gman_ip=10.0.1.2
set gman_port=4019
set nslaves=4
set t_folder=template\%~n0



xcopy /e /q /y %t_folder%\genie.exe host\
xcopy /e /q /y %t_folder% master\
FOR /L %%i IN (1,1,%nslaves%) DO ( xcopy /e /q /y %t_folder% Slave%%i\)

if "%1" == "" ( echo no bootstrap)
else(

 xcopy /e /q /y template\%1 master\
 FOR /L %%i IN (1,1,%nslaves%) DO ( xcopy /e /q /y template\%1 Slave%%i\)
)


REM Genie Host
start "Genie Host:%gman_port%" /d .\host .\host\Genie.exe /port %gman_port% /NFAIL 2 /RUNTIME .2

REM Genie Slaves
FOR /D %%X IN (.\Slave*) DO (
echo %%~nX
 start "Genie %%X: /host %gman_ip%:%gman_port%" /d %%X %%X\Genie.exe /interval 1.0 /console off /name %%~nX /host %gman_ip%:%gman_port%
)

REM UNCOMMENT THE NEXT LINE IF YOU WOULD LIKE TO MONITOR THE.ext FILE
REM start /d .\master fwait.exe storage5.ext external verbose

REM PEST++, etc
start /d .\master .\startPest.bat %~n0 %gman_ip%:%gman_port%