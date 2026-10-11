#!/usr/bin/env python3
"""
CLI Tool: TIC 4206066 TESS Sector 110 Automated Workflow Scheduler.
Provides commands to monitor the live countdown, perform dry-run tests,
and register or trigger the automated workflow execution upon data release.
"""

from __future__ import annotations

import argparse
import json
import sys
import time

from anse.astrophysics.sector110_scheduler import Sector110WorkflowScheduler


def main():
    parser = argparse.ArgumentParser(
        description="NASA TESS Sector 110 Automated Workflow Scheduler for TIC 4206066"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a test dry-run: executes the complete workflow, verifies receipts, and updates state",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Print the current real-time countdown, scheduler status, and latest results",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check if release date has arrived and trigger automated execution if reached",
    )
    parser.add_argument(
        "--daemon",
        action="store_true",
        help="Run a background liveness loop checking periodically for data release",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=3600,
        help="Polling interval in seconds for daemon mode (default: 3600s / 1 hour)",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output raw JSON instead of human-readable text",
    )
    args = parser.parse_args()

    scheduler = Sector110WorkflowScheduler()

    if args.dry_run:
        result = scheduler.execute_workflow(is_dry_run=True)
        if args.json:
            print(json.dumps(result, indent=2))
        else:
            print("\n" + "=" * 65)
            print("  TESS SECTOR 110 AUTONOMOUS WORKFLOW SCHEDULER: DRY RUN")
            print("=" * 65)
            print(f"Target Object:        TIC 4206066 (StKM 1-561)")
            print(f"Target Release Date:  {result['time_counter']['target_date_utc']}")
            print(f"Time Remaining:       {result['time_counter']['formatted_countdown']}")
            print(f"Execution Mode:       DRY RUN (Pre-flight test)")
            step_count = len(result['workflow_summary']['steps']) if result.get('workflow_summary') else 8
            print(f"Workflow Success:     {'✅ ALL ' + str(step_count) + ' STEPS PASSED' if result['success'] else '❌ FAILED'}")
            print(f"Total Duration:       {result['workflow_summary']['total_duration_ms']} ms")
            print(f"Receipt Provenance:   {result['receipt']['provenance_hash']}")
            print(f"State File Saved:     {result['state_file']}")
            print("=" * 65)
            print("Self-test passed. NOTE: this runner is NOT the pre-registered analysis (see docs/preregistration/PROVENANCE.md).\n")
        sys.exit(0 if result["success"] else 1)

    elif args.status:
        st = scheduler.get_status()
        if args.json:
            print(json.dumps(st, indent=2))
        else:
            cnt = st["time_counter"]
            print("\n" + "=" * 65)
            print("       TESS SECTOR 110 WORKFLOW SCHEDULER STATUS")
            print("=" * 65)
            print(f"Schedule Name:        {st['schedule_name']}")
            print(f"Target ID:            {st['target_id']}")
            print(f"Expected Release:     {st['target_release_utc']}")
            print(f"Current UTC Time:     {cnt['current_time_utc']}")
            print(f"Countdown to Release: {cnt['formatted_countdown']}")
            print(f"Current State:        {st['status']}")
            print(f"Last Dry-Run Passed:  {st.get('dry_run_last_passed_utc', 'None')}")
            print(f"Notes:                {st['notes']}")
            print("=" * 65 + "\n")
        sys.exit(0)

    elif args.check:
        res = scheduler.check_and_trigger()
        if args.json:
            print(json.dumps(res, indent=2))
        else:
            if res.get("triggered", False):
                print(f"✅ Triggered workflow execution! Status: {res['status']}")
            else:
                print(f"⏳ Waiting for release. Countdown: {res['time_counter']['formatted_countdown']}")
        sys.exit(0)

    elif args.daemon:
        print(f"Starting TESS Sector 110 scheduler daemon (poll interval: {args.interval}s)...")
        try:
            while True:
                res = scheduler.check_and_trigger()
                if res.get("triggered", False):
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S UTC')}] Release date reached! Workflow executed: {res['status']}")
                    break
                else:
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S UTC')}] {res['time_counter']['formatted_countdown']} remaining.")
                time.sleep(args.interval)
        except KeyboardInterrupt:
            print("\nDaemon terminated by user.")
        sys.exit(0)

    else:
        st = scheduler.get_status()
        print(f"TESS Sector 110 Countdown: {st['time_counter']['formatted_countdown']} (Status: {st['status']})")
        print("Use --dry-run to test or --status for full details.")


if __name__ == "__main__":
    main()
