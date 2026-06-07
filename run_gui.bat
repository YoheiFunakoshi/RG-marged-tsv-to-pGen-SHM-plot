@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\pythonw.exe" (
  echo Python environment not found: .venv
  echo Create it with:
  echo   python -m venv .venv
  echo   .venv\Scripts\python.exe -m pip install -r requirements.txt
  pause
  exit /b 1
)
".venv\Scripts\pythonw.exe" "AIRR_pGen_SHM_plot_app.pyw"
