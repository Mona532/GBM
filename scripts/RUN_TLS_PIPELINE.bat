@echo off
echo ============================================================
echo  SpaLinker TLS Detection Pipeline
echo  Step 1: Extract h5ad data (Python)
echo  Step 2: Run TLS analysis (R)
echo ============================================================
echo.

set PYTHON=C:\Users\qingy\AppData\Local\miniconda3\envs\cell2loc\python.exe
set RSCRIPT=E:\GBM\R\R-4.3.2\bin\Rscript.exe
set R_HOME=E:\GBM\R\R-4.3.2
set R_LIBS_USER=E:\GBM\R\R-4.3.2\library

echo [Step 1/2] Extracting data from h5ad files...
echo   This may take 5-10 minutes for all samples.
%PYTHON% E:\GBM\scripts\01_extract_all_h5ad.py
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: Python extraction failed!
    pause
    exit /b 1
)
echo   Done!
echo.

echo [Step 2/2] Running SpaLinker TLS analysis...
%RSCRIPT% --no-environ E:\GBM\scripts\02_run_spalinker_tls.R
if %ERRORLEVEL% NEQ 0 (
    echo ERROR: R analysis failed!
    pause
    exit /b 1
)
echo   Done!
echo.

echo ============================================================
echo  Pipeline complete!
echo  Results: E:\GBM\results\spalinker_tls_final\
echo ============================================================
pause
