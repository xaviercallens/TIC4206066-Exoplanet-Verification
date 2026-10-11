#!/usr/bin/env python3
"""
TESS Sector 110 Autonomous Workflow Scheduler & Real-Time Countdown Engine.
Orchestrates automated execution of the ExoElenchus-TIC4206066 pipeline when
NASA releases TESS Sector 110 high-cadence data in late November 2026.

Features:
1. Real-time countdown timer (Days, Hours, Minutes, Seconds remaining).
2. Autonomous execution trigger when release timestamp is reached.
3. Full Dry-Run simulation engine to verify end-to-end execution ahead of time.
4. Cryptographic receipt and state persistence in results/tic4206066_scheduler_state.json.
5. Integration with the ANSE Web Command Deck and automated test harness.
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from typing import Any

import numpy as np

from scripts.run_tic4206066_workflow import TIC4206066Workflow
from subprojects.tic4206066_sector110_runner import (
    Sector110PredictionReceipt,
    execute_preregistered_runner,
)

logger = logging.getLogger("Sector110Scheduler")


class NumpyEncoder(json.JSONEncoder):
    """Custom JSON encoder for NumPy types."""
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, (np.bool_, bool)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


# NASA TESS Sector 110 observations: 2026-10-30 to 2026-11-25 (see subprojects runner).
# The public-release timestamp below is a NOMINAL ESTIMATE (observation end + ~3 days),
# not taken from a NASA/MAST schedule; override via target_release_iso when known.
SECTOR_110_DEFAULT_RELEASE_ISO = "2026-11-28T18:00:00+00:00"
STATE_FILE_PATH = "results/tic4206066_scheduler_state.json"


@dataclass
class TimeCounter:
    target_date_utc: str
    current_time_utc: str
    time_remaining_seconds: float
    days: int
    hours: int
    minutes: int
    seconds: int
    is_reached: bool
    formatted_countdown: str


@dataclass
class SchedulerState:
    schedule_name: str
    target_id: str
    target_release_utc: str
    status: str  # ARMED, DRY_RUN_VERIFIED, EXECUTING, COMPLETED, FAILED
    last_checked_utc: str
    time_counter: dict[str, Any]
    latest_workflow_summary: dict[str, Any] | None
    execution_receipt: dict[str, Any] | None
    dry_run_last_passed_utc: str | None
    notes: str


def compute_time_counter(
    target_iso: str = SECTOR_110_DEFAULT_RELEASE_ISO,
    simulated_now: datetime | None = None,
) -> TimeCounter:
    """
    Computes real-time or simulated countdown until NASA TESS Sector 110 data release.
    """
    target_dt = datetime.fromisoformat(target_iso)
    if target_dt.tzinfo is None:
        target_dt = target_dt.replace(tzinfo=UTC)

    now_dt = simulated_now or datetime.now(UTC)
    if now_dt.tzinfo is None:
        now_dt = now_dt.replace(tzinfo=UTC)

    diff = (target_dt - now_dt).total_seconds()
    is_reached = diff <= 0.0

    rem = max(0.0, diff)
    days = int(rem // 86400)
    hours = int((rem % 86400) // 3600)
    mins = int((rem % 3600) // 60)
    secs = int(rem % 60)

    if is_reached:
        formatted = "00d 00h 00m 00s (RELEASE DATE REACHED)"
    else:
        formatted = f"{days:02d}d {hours:02d}h {mins:02d}m {secs:02d}s"

    return TimeCounter(
        target_date_utc=target_dt.isoformat(),
        current_time_utc=now_dt.isoformat(),
        time_remaining_seconds=round(diff, 1),
        days=days,
        hours=hours,
        minutes=mins,
        seconds=secs,
        is_reached=is_reached,
        formatted_countdown=formatted,
    )


class Sector110WorkflowScheduler:
    """Manages scheduling, countdown monitoring, dry runs, and autonomous execution."""

    def __init__(
        self,
        target_release_iso: str = SECTOR_110_DEFAULT_RELEASE_ISO,
        state_path: str = STATE_FILE_PATH,
    ):
        self.target_release_iso = target_release_iso
        self.state_path = state_path
        os.makedirs(os.path.dirname(state_path) or ".", exist_ok=True)
        self.state = self._load_or_init_state()

    def _load_or_init_state(self) -> SchedulerState:
        if os.path.exists(self.state_path):
            try:
                with open(self.state_path) as f:
                    data = json.load(f)
                return SchedulerState(**data)
            except Exception as e:
                logger.warning(f"Could not load state from {self.state_path}: {e}")

        counter = compute_time_counter(self.target_release_iso)
        return SchedulerState(
            schedule_name="NASA TESS Sector 110 Automated Ingestion & Elenchus Vetting",
            target_id="TIC 4206066",
            target_release_utc=self.target_release_iso,
            status="ARMED",
            last_checked_utc=datetime.now(UTC).isoformat(),
            time_counter=asdict(counter),
            latest_workflow_summary=None,
            execution_receipt=None,
            dry_run_last_passed_utc=None,
            notes="Armed for automated execution upon NASA MAST public release in late November 2026.",
        )

    def save_state(self) -> None:
        """Persists scheduler state to disk."""
        with open(self.state_path, "w") as f:
            json.dump(asdict(self.state), f, indent=2, cls=NumpyEncoder)

    def get_status(self) -> dict[str, Any]:
        """Refreshes the live countdown timer and returns current status."""
        counter = compute_time_counter(self.target_release_iso)
        self.state.last_checked_utc = datetime.now(UTC).isoformat()
        self.state.time_counter = asdict(counter)
        self.save_state()
        return asdict(self.state)

    def execute_workflow(self, is_dry_run: bool = False) -> dict[str, Any]:
        """
        Executes the full ExoElenchus-TIC4206066 workflow and updates receipts and state.
        """
        logger.info(f"Executing workflow (dry_run={is_dry_run})...")
        self.state.status = "EXECUTING"
        self.save_state()

        # Step A: Run complete 7-step reproducible workflow
        workflow = TIC4206066Workflow()
        workflow_summary = workflow.run_all()

        # Step B: Run and sign Sector 110 Pre-Registration receipt
        receipt: Sector110PredictionReceipt = execute_preregistered_runner()

        now_str = datetime.now(UTC).isoformat()
        self.state.latest_workflow_summary = workflow_summary
        self.state.execution_receipt = asdict(receipt)
        self.state.last_checked_utc = now_str

        if is_dry_run:
            self.state.status = "DRY_RUN_VERIFIED"
            self.state.dry_run_last_passed_utc = now_str
            self.state.notes = (
                "Dry run: pipeline executed on archived Sectors 6/32/98 and a SYNTHETIC "
                "Sector 110 self-test. No Sector 110 data were evaluated. "
                f"workflow_all_passed={bool(workflow_summary.get('all_passed'))}."
            )
            live_outcome: dict[str, Any] | None = None
        else:
            live_outcome = self._evaluate_live_sector_110()
            if live_outcome is None:
                self.state.status = "AWAITING_DATA"
                self.state.notes = (
                    f"{now_str}: release date reached but no Sector 110 light-curve FITS was "
                    "found in data/astrophysics/tic4206066/mast. Nothing was evaluated; this is "
                    "NOT a verification result."
                )
            else:
                self.state.status = "COMPLETED" if workflow_summary.get("all_passed") else "FAILED"
                self.state.notes = (
                    f"{now_str}: Sector 110 evaluated on real data. "
                    f"Signal A outcome: {live_outcome['signal_a_outcome']}."
                )

        counter = compute_time_counter(self.target_release_iso)
        self.state.time_counter = asdict(counter)
        self.save_state()

        success = bool(workflow_summary.get("all_passed", False)) and (
            is_dry_run or live_outcome is not None
        )
        return {
            "success": success,
            "status": self.state.status,
            "is_dry_run": is_dry_run,
            "time_counter": asdict(counter),
            "workflow_summary": workflow_summary,
            "receipt": asdict(receipt),
            "live_sector_110_outcome": live_outcome,
            "state_file": self.state_path,
        }

    @staticmethod
    def _evaluate_live_sector_110() -> dict[str, Any] | None:
        """Evaluate real Sector 110 photometry if (and only if) a FITS file is present."""
        from subprojects.tic4206066_sector110_runner import (
            evaluate_sector_110_data,
            find_sector110_fits,
            load_sector110_fits,
        )

        path = find_sector110_fits()
        if path is None:
            return None
        lc = load_sector110_fits(path)
        out = evaluate_sector_110_data(lc["times"], lc["fluxes"], lc["pixel_x"], lc["pixel_y"])
        out["source_file"] = path
        return out

    def check_and_trigger(self) -> dict[str, Any]:
        """
        Polls the countdown timer. If target release timestamp has arrived,
        automatically triggers the live workflow execution.
        """
        counter = compute_time_counter(self.target_release_iso)
        self.state.time_counter = asdict(counter)
        self.state.last_checked_utc = datetime.now(UTC).isoformat()

        if counter.is_reached:
            logger.info("Target release timestamp reached! Launching automated workflow...")
            return self.execute_workflow(is_dry_run=False)

        self.save_state()
        return {
            "triggered": False,
            "status": self.state.status,
            "time_counter": asdict(counter),
            "message": f"Waiting for release date {self.target_release_iso}. {counter.formatted_countdown}",
        }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="TESS Sector 110 Workflow Scheduler")
    parser.add_argument("--dry-run", action="store_true", help="Execute test dry-run")
    parser.add_argument("--status", action="store_true", help="Get current countdown and state")
    parser.add_argument("--check", action="store_true", help="Check date and trigger if reached")
    args = parser.parse_args()

    scheduler = Sector110WorkflowScheduler()

    if args.dry_run:
        res = scheduler.execute_workflow(is_dry_run=True)
        print("\n=== TESS SECTOR 110 WORKFLOW SCHEDULER: DRY RUN ===")
        print(f"Status: {res['status']}")
        print(f"Countdown: {res['time_counter']['formatted_countdown']}")
        print(f"Workflow Passed: {res['success']}")
        print(f"Total Workflow Duration: {res['workflow_summary']['total_duration_ms']} ms")
        print(f"Execution Receipt SHA-256: {res['receipt']['provenance_hash'][:16]}...")
        print(f"State saved to {res['state_file']}\n")
    elif args.status:
        st = scheduler.get_status()
        print("\n=== TESS SECTOR 110 SCHEDULER STATUS ===")
        print(f"Target: {st['target_id']} (Release: {st['target_release_utc']})")
        print(f"Countdown: {st['time_counter']['formatted_countdown']}")
        print(f"Scheduler Status: {st['status']}")
        print(f"Notes: {st['notes']}\n")
    elif args.check:
        res = scheduler.check_and_trigger()
        print(json.dumps(res, indent=2))
    else:
        st = scheduler.get_status()
        print(f"Countdown until Sector 110 data release: {st['time_counter']['formatted_countdown']}")
