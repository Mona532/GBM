@echo off
setlocal

set RSCRIPT=E:\GBM\R\R-4.3.2\bin\Rscript.exe

echo [1/3] Build TLS component pseudobulk
"%RSCRIPT%" "E:\GBM\scripts\03_build_tls_component_pseudobulk.R"
if errorlevel 1 goto :fail

echo [2/3] Run TLS component NMF
"%RSCRIPT%" "E:\GBM\scripts\04_run_tls_component_nmf.R"
if errorlevel 1 goto :fail

echo [3/3] Run TLS component markers
"%RSCRIPT%" "E:\GBM\scripts\05b_run_tls_component_markers_aligned.R"
if errorlevel 1 goto :fail

echo TLS component ecotype pipeline completed.
exit /b 0

:fail
echo Pipeline failed with exit code %errorlevel%.
exit /b %errorlevel%
