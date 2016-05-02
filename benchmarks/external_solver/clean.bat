@echo off
REM CLEAN UP
rmdir /S /Q .\host
rmdir /S /Q .\master
FOR /D %%X IN (.\slave*) DO (rmdir /S /Q %%X)