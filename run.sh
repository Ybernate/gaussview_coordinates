#!/usr/bin/env bash
# Convenience launcher for MolGauss.
# Creates the virtualenv and installs dependencies on first run, then starts
# the Streamlit app.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d ".venv" ]; then
  echo "Creating virtual environment…"
  python3 -m venv .venv
  # shellcheck disable=SC1091
  source .venv/bin/activate
  pip install --upgrade pip
  pip install -r requirements.txt
else
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

exec streamlit run app.py "$@"
