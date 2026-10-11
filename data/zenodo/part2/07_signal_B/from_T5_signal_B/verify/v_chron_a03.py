"""File-time chronology including the frozen A03 harmonic audit (read-only os.stat), UTC."""
import os, json, datetime as dt
ROOT = "/Users/paulsmith/space <3/tess-jev/"
fs = ["config/amendments_registry.jsonl", "data/harmonic/qlp_downloads.jsonl", "logs/harmonic_run.log",
      "data/harmonic/per_candidate/tic4206066_p11.132976.json", "data/harmonic/harmonic_audit.csv",
      "candidates/tic4206066_p11.132976/harmonic.json", "addenda/pool_multisector/multisector_results.json",
      "addenda/pool_multisector/SUMMARY.md", "addenda/deep_followup_tic4206066/signal_11d.py",
      "addenda/deep_followup_tic4206066/out/signal_11d.json"]
u = lambda s: dt.datetime.fromtimestamp(s, dt.timezone.utc).isoformat()
out = {f: dict(birth_utc=u(os.stat(ROOT + f).st_birthtime), mtime_utc=u(os.stat(ROOT + f).st_mtime)) for f in fs}
h = json.load(open(ROOT + "data/harmonic/per_candidate/tic4206066_p11.132976.json"))
out["A03_B_events_at_P"] = [dict(sector=e["sector"], k=e["k"], t_mid=round(e["t_mid"], 3), cls=e["class"], snr=e["snr"], present=e["present"], contradicting=e["contradicting"]) for e in h["hypotheses"]["P"]["events"]]
out["A03_B_status"] = dict(period_status=h["period_status"], supported=h.get("supported_hypotheses"), n_contra_P=h.get("n_contradicting_epochs_at_P"), sector_shifts_d=h["hypotheses"]["P"]["sector_shifts_d"])
json.dump(out, open(ROOT + "addenda/tic4206066_FINAL_RECONCILIATION/_work/T5_signal_B/verify/v_chron_a03.json", "w"), indent=1)
for k, v in out.items(): print(k, v)
