@echo off
setlocal enabledelayedexpansion
title DNA Sequence Analyzer - Build Standalone Executable
color 0A

echo ======================================================================
echo           DNA Sequence Analyzer - Executable Builder
echo ======================================================================
echo.
echo Checking Python environment...

REM Check for Python launcher or python command
set PYTHON_CMD=
py --version >nul 2>&1
if %errorlevel% equ 0 (
    set PYTHON_CMD=py
    goto :PYTHON_FOUND
)

python --version >nul 2>&1
if %errorlevel% equ 0 (
    set PYTHON_CMD=python
    goto :PYTHON_FOUND
)

color 0C
echo [ERROR] Python was not found on your system.
echo Please install Python 3.9 or newer from https://www.python.org/
echo Make sure to check "Add Python to PATH" during installation.
echo.
pause
exit /b 1

:PYTHON_FOUND
echo [OK] Using Python command: %PYTHON_CMD%
%PYTHON_CMD% --version

echo.
echo Checking and installing build dependencies (matplotlib, pyinstaller)...
%PYTHON_CMD% -m pip install --quiet --upgrade pip
%PYTHON_CMD% -m pip install --quiet matplotlib pyinstaller
if %errorlevel% neq 0 (
    color 0C
    echo [ERROR] Failed to install required build dependencies.
    pause
    exit /b 1
)
echo [OK] Dependencies ready.

echo.
echo Running automated unit tests before building...
%PYTHON_CMD% -m unittest discover tests
if %errorlevel% neq 0 (
    color 0C
    echo [ERROR] Unit tests failed! Aborting executable build.
    pause
    exit /b 1
)
echo [OK] All automated tests passed.

echo.
echo Bundling standalone application with PyInstaller...
%PYTHON_CMD% build_exe.py
if %errorlevel% neq 0 (
    color 0C
    echo [ERROR] Packaging failed.
    pause
    exit /b 1
)

color 0A
echo.
echo ======================================================================
echo [SUCCESS] Standalone executable ready!
echo File: DNA Sequence Analyzer.exe
echo Double-click "DNA Sequence Analyzer.exe" to launch the application.
echo ======================================================================
echo.
pause
