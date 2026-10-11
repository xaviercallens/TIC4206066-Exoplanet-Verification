#!/bin/sh
# Rebuild everything in this package. Uses the project's .venv read-only.
set -e
cd "$(dirname "$0")"
PY="../../../.venv/bin/python"
"$PY" s01_inputs.py
"$PY" s02_before_global_detrend.py
"$PY" s03_after_decorrelated.py
"$PY" s03b_t9_method_replay.py
"$PY" s04_b_s32_events.py
"$PY" s05_plots.py
"$PY" s06_hashes.py
