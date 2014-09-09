@echo off
setlocal EnableDelayedExpansion
set count=0
FOR /F "delims=" %%a in (%1) do (
    set /A count+= 1
    set var!count!=%%a
)
set "var=!var:*:=!"

cd !var1!
"C:\RZWQM2\bin\pest\par2par" .\pest\par2par.in
"C:\RZWQM2\bin\rzwqmrelease.exe"
cd !var2!
endlocal
@echo on