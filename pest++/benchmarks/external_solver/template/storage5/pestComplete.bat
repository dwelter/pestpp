@echo off
set RETURN=END 
set ARCH_DEL=NOOP
REM set ARCH_DEL=DELETE
REM set ARCH_DEL=ARCHIVE

if exist %1.ext (

 if exist %1 (
  set RETURN=CONTINUECOMPLETE
  REM GOTO DELETE
  GOTO %ARCH_DEL%
 )

:CONTINUECOMPLETE
 if NOT exist %1 ( mkdir %1 )

 pbin2ascii.exe %1.ext %1.runs %1/RUN
 echo runs_are_needed > %1.ctl
 start "fwait: %1.ctl" /d . .\fwait.exe %1.ctl runs_are_done "restartPest.bat %1" verbose
 GOTO END

) else (

 set RETURN=END
 echo all_done > %1.ctl
 GOTO %ARCH_DEL%
)

:NOOP
REM DO NOTHING
GOTO %RETURN%

:DELETE
REM TO REMOVE THE FILES FOR PAST RUNS...
del /q /f %1\*.run
del /q /f %1\*.par
GOTO %RETURN%



:ARCHIVE
REM TO SAVE THEM SOMEWHERE..
setlocal enableDelayedExpansion

pushd .
cd %1

set "baseName=processed"
set "n=0"
for /f "delims=" %%F in (
  '2^>nul dir /b /ad "%baseName%*."^|findstr /xri "%baseName%[0-9]*"'
) do (
  set "name=%%F"
  set "name=!name:*%baseName%=!"
  if !name! gtr !n! set "n=!name!"
)
set /a n+=1

popd

mkdir %1\%baseName%%n%

copy /a /y %1\*.run %1\%baseName%%n%
copy /a /y %1\*.par %1\%baseName%%n%
copy /a /y %1.runs %1\%baseName%%n%
copy /a /y %1.results %1\%baseName%%n%

GOTO DELETE



:END
exit
