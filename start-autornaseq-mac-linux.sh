#!/usr/bin/env sh
set -e
cd "$(dirname "$0")"
docker compose up -d --build shinyapp
printf '\nAutoRNAseq is starting.\n'
printf 'Open http://localhost:3838 in your browser.\n'
printf 'To stop it later, run: docker compose down\n'
