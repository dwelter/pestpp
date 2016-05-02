@echo off
ascii2pbin.exe %1.ext %1.results
del /q /f %1.ext
pest++.exe %1 /e /r

start /d . .\pestComplete.bat %1