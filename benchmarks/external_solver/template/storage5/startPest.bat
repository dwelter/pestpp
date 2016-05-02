@echo off

rmdir /S /Q %1
del /q /f %1.runs
del /q /f %1.results
del /q /f %1.ctl
del /q /f %1.ext

start "p-sensan: /h %2" /d . .\p-sensan.exe /f %1.scf /h %2 /e %1 /ef %1

pest++.exe %1 /e

start /d . .\pestComplete.bat %1
