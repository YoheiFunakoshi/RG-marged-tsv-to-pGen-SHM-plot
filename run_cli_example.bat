@echo off
cd /d "%~dp0"
if "%~1"=="" (
  echo Usage: run_cli_example.bat path\to\sample.igblast.airr.tsv
  pause
  exit /b 1
)
".venv\Scripts\python.exe" "airr_pgen_shm_plot_beta1.py" --input "%~1" --outdir "%~dp1" --sample "%~n1"
pause
