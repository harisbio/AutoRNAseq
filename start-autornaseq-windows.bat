@echo off
setlocal
cd /d %~dp0
docker compose up -d --build shinyapp
echo.
echo AutoRNAseq is starting.
echo Open http://localhost:3838 in your browser.
echo To stop it later, run: docker compose down
endlocal
