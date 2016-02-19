@echo off
REM This part of batch file added by SVDAPREP
REM
REM Delete model input files.
del del C:\PEST\slave1\par2par_rzdat.in > nul
REM
REM Run PARCALC to compute base parameters from super parameters.
C:\PEST\bin\parcalc > nul
REM
REM The following is copied directly from file C:\PEST\slave1\rzwqmb.bat
REM
cd C:\RZWQM2\Nashu_Fox2\nashua13
C:\PEST\bin\par2par   C:\PEST\slave1\par2par_rzdat.in

c:\rzwqm2\bin\rzwqmrelease.exe

