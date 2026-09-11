@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if defined BSIM_ROOT (
  set "CP=%BSIM_ROOT%\dist\build;%BSIM_ROOT%\lib\core.jar;%BSIM_ROOT%\lib\vecmath.jar;%BSIM_ROOT%\lib\objimport.jar"
) else (
  set "CP=..\..\dist\build;..\..\lib\core.jar;..\..\lib\vecmath.jar;..\..\lib\objimport.jar"
)

echo Compiling PocketNeckT3...
javac -encoding UTF-8 -cp "%CP%" -d . bsim\BSimPocketNeckT3.java
if errorlevel 1 exit /b 1

if "%~1"=="" (
  echo.
  echo Compiled. From this folder run for example:
  echo   compile_and_run.cmd config\field_w100.properties
  echo   compile_and_run.cmd config\live_w20.properties
  echo   compile_and_run.cmd config\live_w20.properties preview
  exit /b 0
)

echo Running PocketNeckT3.BSimPocketNeckT3 %*
java -cp ".;%CP%" PocketNeckT3.BSimPocketNeckT3 %*
exit /b %ERRORLEVEL%
