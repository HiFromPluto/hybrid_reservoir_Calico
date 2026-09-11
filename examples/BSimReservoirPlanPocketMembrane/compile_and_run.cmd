@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if defined BSIM_ROOT (
  set "CP=%BSIM_ROOT%\dist\build;%BSIM_ROOT%\lib\core.jar;%BSIM_ROOT%\lib\vecmath.jar;%BSIM_ROOT%\lib\objimport.jar"
) else (
  set "CP=..\..\dist\build;..\..\lib\core.jar;..\..\lib\vecmath.jar;..\..\lib\objimport.jar"
)

echo Compiling PocketMembrane...
javac -encoding UTF-8 -cp "%CP%" -d . bsim\BSimPocketMembrane.java
if errorlevel 1 exit /b 1

if "%~1"=="" (
  echo.
  echo Compiled. From this folder run for example:
  echo   compile_and_run.cmd config\live_w20_mem.properties
  echo   compile_and_run.cmd config\live_w20_mem.properties preview
  exit /b 0
)

echo Running PocketMembrane.BSimPocketMembrane %*
java -cp ".;%CP%" PocketMembrane.BSimPocketMembrane %*
exit /b %ERRORLEVEL%
