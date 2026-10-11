"""
ExoElenchus: TIC 4206066 Discovery & Verification Portal.
Cloud Run ASGI Service serving the Stitch-designed astrophysics dashboard,
live TESS Sector 110 countdown ticker, and autonomous verification pipeline.
"""

from __future__ import annotations

from datetime import datetime, timezone
import json
import logging
import math
import os
from pathlib import Path
import time
from typing import Any
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
import numpy as np

logger = logging.getLogger("exoelenchus.portal")

app = FastAPI(
    title="ExoElenchus TIC 4206066 Discovery Portal",
    description="Astrophysics portal for candidate exoplanet TIC 4206066 b and TESS Sector 110 verification.",
    version="1.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
FORMAL_DIR = BASE_DIR / "formal"
TEMPLATE_FILE = BASE_DIR / "stitch_template.html"

# Nominal estimate (observation end + ~3 d); NOT an official MAST release schedule.
SECTOR_110_RELEASE_ISO = "2026-11-28T18:00:00+00:00"


def compute_countdown(target_iso: str = SECTOR_110_RELEASE_ISO) -> dict[str, Any]:
    target_dt = datetime.fromisoformat(target_iso)
    now_dt = datetime.now(timezone.utc)
    diff = (target_dt - now_dt).total_seconds()
    is_reached = diff <= 0

    rem = max(0.0, diff)
    days = int(rem // 86400)
    hours = int((rem % 86400) // 3600)
    mins = int((rem % 3600) // 60)
    secs = int(rem % 60)

    formatted = "00d 00h 00m 00s (RELEASE REACHED)" if is_reached else f"{days:02d}d {hours:02d}h {mins:02d}m {secs:02d}s"

    return {
        "target_date_utc": target_dt.isoformat(),
        "current_time_utc": now_dt.isoformat(),
        "time_remaining_seconds": round(diff, 1),
        "days": days,
        "hours": hours,
        "minutes": mins,
        "seconds": secs,
        "is_reached": is_reached,
        "formatted_countdown": formatted,
    }


def _load_json(filename: str | Path) -> dict[str, Any] | None:
    if isinstance(filename, Path) and filename.is_absolute():
        paths = [filename]
    else:
        name = Path(filename).name
        paths = [
            RESULTS_DIR / name,
            DATA_DIR / name,
            DATA_DIR / "astrophysics" / "tic4206066" / name,
            BASE_DIR / name,
        ]
    for p in paths:
        if p.exists():
            try:
                with open(p, encoding="utf-8") as f:
                    return json.load(f)
            except (OSError, ValueError) as exc:
                logger.warning("Failed to parse %s: %s", p, exc)
    return None


def _lean_status(summary: dict[str, Any]) -> dict[str, Any]:
    """Derive the Lean 4 status from the recorded workflow step 8 (never hardcoded)."""
    step = next((s for s in summary.get("steps", []) if s.get("step_id") == 8), None)
    if step is None:
        return {"status": "UNKNOWN", "detail": "step 8 absent from workflow summary"}
    d = step.get("details", {})
    return {
        "status": "ELABORATED_ZERO_SORRY" if step.get("passed") else "FAILED",
        "zero_sorry": d.get("zero_sorry_soundness"),
        "elaboration_clean": d.get("lean_elaboration_clean"),
        "theorems": d.get("formalized_theorems", []),
        "scope": "Lean checks the stated algebraic invariants only; it does not validate the data.",
    }


@app.get("/", response_class=HTMLResponse)
def get_dashboard() -> str:
    """Serve the Stitch-designed discovery & verification dashboard."""
    if TEMPLATE_FILE.exists():
        return TEMPLATE_FILE.read_text(encoding="utf-8")
    return "<h1>ExoElenchus TIC 4206066 Discovery Portal</h1>"


@app.get("/health")
def health() -> dict[str, str]:
    return {"status": "healthy", "service": "exoplanet-tic4206066-portal"}


@app.get("/api/astrophysics/tic4206066/scheduler")
def get_scheduler_status() -> dict[str, Any]:
    """Returns the live countdown and armed scheduler status."""
    countdown = compute_countdown()
    data = _load_json("tic4206066_scheduler_state.json")
    if data:
        data["time_counter"] = countdown
        data["last_checked_utc"] = datetime.now(timezone.utc).isoformat()
        return data

    return {
        "schedule_name": "NASA TESS Sector 110 Automated Ingestion & Elenchus Vetting",
        "target_id": "TIC 4206066",
        "target_release_utc": SECTOR_110_RELEASE_ISO,
        "status": "STATE_UNAVAILABLE",
        "time_counter": countdown,
        "notes": "No scheduler state file found; no verification status is claimed. "
        "Release date is a nominal estimate, not an official MAST schedule.",
    }


@app.post("/api/astrophysics/tic4206066/scheduler/dry-run")
@app.get("/api/astrophysics/tic4206066/verify")
def execute_dry_run() -> JSONResponse:
    """
    Return the RECORDED output of the verification workflow (nothing is re-run or invented).
    If the recorded artifacts are missing the endpoint answers 503 instead of fabricating a result.
    """
    summary = _load_json("tic4206066_workflow_summary.json")
    receipt = _load_json("tic4206066_execution_receipt.json")
    if summary is None or receipt is None:
        return JSONResponse(
            status_code=503,
            content={
                "success": False,
                "status": "ARTIFACTS_UNAVAILABLE",
                "detail": "workflow summary / execution receipt not found; run scripts/run_tic4206066_workflow.py",
            },
        )

    return JSONResponse(
        content={
            "success": bool(summary.get("all_passed")),
            "status": receipt.get("overall_status"),
            "is_dry_run": True,
            "data_kind": "archived Sectors 6/32/98 + synthetic Sector 110 self-test; no Sector 110 data evaluated",
            "time_counter": compute_countdown(),
            "workflow_summary": summary,
            "receipt": receipt,
            "provenance_hash": receipt.get("provenance_hash"),
            "lean4_proof": _lean_status(summary),
            "git_hash": summary.get("git_hash"),
            "seed": summary.get("seed"),
            "references": {
                "zenodo_preprint_doi": "10.5281/zenodo.22967456",
                "zenodo_data_code_doi": "10.5281/zenodo.22967099",
                "zenodo_preregistration_doi": "10.5281/zenodo.23175179",
                "nasa_mission": "TESS (Transiting Exoplanet Survey Satellite)",
                "archive": "MAST STScI",
            },
        }
    )


@app.get("/api/astrophysics/tic4206066/status")
def get_status() -> dict[str, Any]:
    """Retrieve certified execution receipt and metadata for TIC 4206066."""
    receipt = _load_json("tic4206066_execution_receipt.json")
    meta = _load_json("metadata.json")
    if receipt:
        result = dict(receipt)
        if meta:
            result["metadata"] = meta
        return result
    if meta:
        return {"target_id": "TIC 4206066", "metadata": meta, "status": "METADATA_ONLY"}
    return {"target_id": "TIC 4206066", "status": "UNKNOWN"}


@app.get("/api/astrophysics/tic4206066/lean")
def get_lean_proof_details() -> dict[str, Any]:
    """Returns the Lean 4 formal proof verification details and theorem statements."""
    summary = _load_json("tic4206066_workflow_summary.json") or {}
    lean_file = FORMAL_DIR / "ANSE" / "ExoplanetTransit.lean"
    source = ""
    if lean_file.exists():
        try:
            source = lean_file.read_text(encoding="utf-8")
        except OSError:
            pass

    return {
        "status": _lean_status(summary),
        "lean_toolchain": "leanprover/lean4:v4.34.0-rc2",
        "formal_source_available": bool(source),
        "theorem_count": len(_lean_status(summary).get("theorems", [])),
    }


@app.post("/api/astrophysics/tic4206066/reproduce")
def reproduce_tic4206066() -> dict[str, Any]:
    """
    Executes or returns the Socratic Elenchus verification results on real archived data.
    """
    try:
        from anse.astrophysics.ephemeris import load_ephemeris_a
        from anse.astrophysics.real_data_audit import verify_consumed_zenodo_files
        from anse.astrophysics.real_data_ingest import audit_real_data_sources

        eph = load_ephemeris_a()
        real_audit = audit_real_data_sources()
        zenodo_check = verify_consumed_zenodo_files()

        return {
            "status": "success",
            "target": "TIC 4206066",
            "candidate": "TIC 4206066.01",
            "period_days": eph.period_d,
            "ephemeris_source": eph.source,
            "data_integrity": zenodo_check,
            "real_audit": real_audit,
            "time_counter": compute_countdown(),
        }
    except Exception as exc:
        logger.error("Reproduction error: %s", exc)
        return {
            "status": "error",
            "detail": str(exc),
            "fallback": _load_json("tic4206066_workflow_summary.json"),
        }


if __name__ == "__main__":
    import uvicorn
    port = int(os.environ.get("PORT", 8080))
    uvicorn.run("app:app", host="0.0.0.0", port=port, reload=False)
