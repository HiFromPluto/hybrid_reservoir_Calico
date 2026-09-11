@echo off
setlocal EnableExtensions
cd /d "%~dp0"

if defined BSIM_ROOT (
  set "ROOT=%BSIM_ROOT%"
) else (
  set "ROOT=%~dp0..\.."
)

set "CP=%ROOT%\dist\build;%ROOT%\lib\core.jar;%ROOT%\lib\vecmath.jar;%ROOT%\lib\objimport.jar"

if not exist "%ROOT%\dist\build\bsim\BSim.class" (
  echo Compiling BSim library into dist\build...
  if not exist "%ROOT%\dist\build" mkdir "%ROOT%\dist\build"
  dir /s /b "%ROOT%\src\bsim\*.java" | findstr /v /i "\\tests\\" > "%TEMP%\bsim_job2_sources.txt"
  javac -encoding UTF-8 -cp "%ROOT%\lib\core.jar;%ROOT%\lib\vecmath.jar;%ROOT%\lib\objimport.jar" -d "%ROOT%\dist\build" @"%TEMP%\bsim_job2_sources.txt"
  if errorlevel 1 exit /b 1
)

echo Compiling BacteriumFromScratch...
javac -encoding UTF-8 -cp "%CP%" -d . ChassisParameters.java EcoliRodCell.java ValdezHertzian.java NutrientField.java Job3Sims.java Job3bSims.java BacteriumFromScratch.java Job5Sims.java SpatialAhlLane.java Job5bSims.java StarvationViability.java Job6Sims.java RunTumbleMotility.java Job7Sims.java
if errorlevel 1 exit /b 1

if /i "%~1"=="job5" (
  echo Running BacteriumFromScratch.Job5Sims
  java -cp ".;%CP%" BacteriumFromScratch.Job5Sims
  exit /b %ERRORLEVEL%
)

if /i "%~1"=="job5b" (
  echo Running BacteriumFromScratch.Job5bSims
  java -cp ".;%CP%" BacteriumFromScratch.Job5bSims
  exit /b %ERRORLEVEL%
)

if /i "%~1"=="job6" (
  echo Running BacteriumFromScratch.Job6Sims
  java -cp ".;%CP%" BacteriumFromScratch.Job6Sims
  exit /b %ERRORLEVEL%
)

if /i "%~1"=="job7" (
  echo Running BacteriumFromScratch.Job7Sims
  java -cp ".;%CP%" BacteriumFromScratch.Job7Sims
  exit /b %ERRORLEVEL%
)

echo Running BacteriumFromScratch.BacteriumFromScratch %*
java -cp ".;%CP%" BacteriumFromScratch.BacteriumFromScratch %*
exit /b %ERRORLEVEL%
