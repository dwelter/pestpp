for /f "tokens=1 delims=," %%a in ('tasklist /fo csv ^|FINDSTR /I /C:"genie.exe"') do call :killprocess %%a
for /f "tokens=1 delims=," %%a in ('tasklist /fo csv ^|FINDSTR /I /C:"pest++.exe"') do call :killprocess %%a
for /f "tokens=1 delims=," %%a in ('tasklist /fo csv ^|FINDSTR /I /C:"p-sensan.exe"') do call :killprocess %%a
for /f "tokens=1 delims=," %%a in ('tasklist /fo csv ^|FINDSTR /I /C:"fwait.exe"') do call :killprocess %%a
for /f "tokens=1 delims=," %%a in ('tasklist /fo csv ^|FINDSTR /I /C:"hydro.exe"') do call :killprocess %%a

:killprocess
echo. |set /p d=killing %*...
taskkill /f /im "%*">nul 2>&1
set err=%errorlevel%
set success=Success
if not %err%==0 set success=fail (err code: %err%)
if %err%==128 set success=fail (process not found)
echo %success%&goto :eof