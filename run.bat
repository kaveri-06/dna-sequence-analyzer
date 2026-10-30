@echo off
setlocal enabledelayedexpansion
title DNA Sequence Analyzer Launcher

REM 1. If standalone executable exists, launch it directly!
if exist "DNA Sequence Analyzer.exe" (
    echo Launching DNA Sequence Analyzer...
    start "" "DNA Sequence Analyzer.exe"
    exit /b 0
)

if exist "dist\DNA Sequence Analyzer.exe" (
    echo Launching DNA Sequence Analyzer...
    start "" "dist\DNA Sequence Analyzer.exe"
    exit /b 0
)

REM 2. Otherwise run using Python
echo ======================================================================
echo                DNA Sequence Analyzer Launcher
echo ======================================================================
echo Standalone executable not found. Starting with Python...
echo.

REM Locate Python
set PYTHON_CMD=
py --version >nul 2>&1
if %errorlevel% equ 0 (
    set PYTHON_CMD=py
    goto :PYTHON_READY
)

python --version >nul 2>&1
if %errorlevel% equ 0 (
    set PYTHON_CMD=python
    goto :PYTHON_READY
)

echo [ERROR] Python is not detected on your computer.
echo.
echo To run this application:
echo 1. Download and install Python from https://www.python.org/
echo 2. During installation, CHECK the box: "Add Python to PATH"
echo 3. Alternatively, run build.bat on a computer with Python to produce the standalone .exe
echo.
pause
exit /b 1

:PYTHON_READY
REM Check if virtual environment exists; if not, create one
if not exist ".venv\Scripts\python.exe" (
    echo Setting up dedicated Python virtual environment (.venv)...
    %PYTHON_CMD% -m venv .venv
    if %errorlevel% neq 0 (
        echo [WARNING] Could not create virtual environment. Running with global Python.
        goto :INSTALL_DEPS_GLOBAL
    )
)

set VENV_PY=.venv\Scripts\python.exe
echo Checking application dependencies...
%VENV_PY% -m pip install --quiet matplotlib
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install required dependencies into virtual environment.
    pause
    exit /b 1
)

echo Launching DNA Sequence Analyzer...
start "" "%VENV_PY%" app.py
exit /b 0

:INSTALL_DEPS_GLOBAL
echo Checking application dependencies...
%PYTHON_CMD% -m pip install --quiet matplotlib
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install matplotlib.
    pause
    exit /b 1
)
echo Launching DNA Sequence Analyzer...
start "" "%PYTHON_CMD%" app.py
exit /b 0
