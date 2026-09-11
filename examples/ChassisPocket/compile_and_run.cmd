@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if defined BSIM_ROOT (
  set "ROOT=%BSIM_ROOT%"
) else (
  set "ROOT=%~dp0..\.."
)

set "CHASSIS=%ROOT%\examples\BacteriumFromScratch"
set "CP=%ROOT%\dist\build;%ROOT%\lib\core.jar;%ROOT%\lib\vecmath.jar;%ROOT%\lib\objimport.jar"

if not exist "%ROOT%\dist\build\bsim\BSim.class" (
  echo Compiling BSim library into dist\build...
  if not exist "%ROOT%\dist\build" mkdir "%ROOT%\dist\build"
  dir /s /b "%ROOT%\src\bsim\*.java" | findstr /v /i "\\tests\\" > "%TEMP%\bsim_i0_sources.txt"
  javac -encoding UTF-8 -cp "%ROOT%\lib\core.jar;%ROOT%\lib\vecmath.jar;%ROOT%\lib\objimport.jar" -d "%ROOT%\dist\build" @"%TEMP%\bsim_i0_sources.txt"
  if errorlevel 1 exit /b 1
)

echo Compiling BacteriumFromScratch chassis (import, not copy)...
javac -encoding UTF-8 -cp "%CP%" -d "%CHASSIS%" "%CHASSIS%\ChassisParameters.java" "%CHASSIS%\EcoliRodCell.java" "%CHASSIS%\ValdezHertzian.java"
if errorlevel 1 exit /b 1

if /I "%~1"=="i0c" (
  echo Compiling ChassisPocketI0c...
  javac -encoding UTF-8 -cp "%CP%;%CHASSIS%" -d . ChassisPocketI0c.java
  if errorlevel 1 exit /b 1
  echo Running ChassisPocket.ChassisPocketI0c %~2
  java -cp ".;%CP%;%CHASSIS%" ChassisPocket.ChassisPocketI0c %~2
  exit /b %ERRORLEVEL%
)

if /I "%~1"=="i0b" (
  echo Compiling ChassisPocketI0b...
  javac -encoding UTF-8 -cp "%CP%;%CHASSIS%" -d . ChassisPocketI0b.java
  if errorlevel 1 exit /b 1
  echo Running ChassisPocket.ChassisPocketI0b %~2
  java -cp ".;%CP%;%CHASSIS%" ChassisPocket.ChassisPocketI0b %~2
  exit /b %ERRORLEVEL%
)

echo Compiling ChassisPocketI0...
javac -encoding UTF-8 -cp "%CP%;%CHASSIS%" -d . ChassisPocketI0.java
if errorlevel 1 exit /b 1

echo Running ChassisPocket.ChassisPocketI0 %*
java -cp ".;%CP%;%CHASSIS%" ChassisPocket.ChassisPocketI0 %*
exit /b %ERRORLEVEL%
