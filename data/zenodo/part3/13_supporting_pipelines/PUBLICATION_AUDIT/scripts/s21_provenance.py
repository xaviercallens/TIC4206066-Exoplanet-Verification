"""Step 21: DATA_PROVENANCE.csv - one row per input file used anywhere in the audit (plus the earlier-session files that
were re-used without re-download, marked as such). Built from data/provenance.jsonl, the FITS headers and the checks of
s01/s02/s06."""
import csv
import json
import os
import sys

from astropy.io import fits

sys.path.insert(0, os.path.dirname(__file__))
from audit_common import AUD, OLD_COMP, OLD_DEEP, sha256_file  # noqa

prov = [json.loads(l) for l in open(os.path.join(AUD, "data", "provenance.jsonl"))]
last = {}
for p in prov:
    if "sha256" in p:
        last[p["dest"]] = p
dl = {r["productFilename"]: r for r in json.load(open(os.path.join(AUD, "data/mast/download_check.json")))}
tc = {c["file"]: c for c in json.load(open(os.path.join(AUD, "data/tesscut/tesscut_check.json")))["cutouts"]}
rows = []


def fits_meta(path):
    try:
        h = fits.open(path)
        hd = {**h[0].header, **h[1].header}
        return dict(release=";".join(str(x) for x in (hd.get("DATA_REL") and f"DATA_REL={hd.get('DATA_REL')}", hd.get("PROCVER") and f"PROCVER={hd.get('PROCVER')}",
                                                       hd.get("CALIB") and f"CALIB={hd.get('CALIB')}", hd.get("TICVER") and f"TICVER={hd.get('TICVER')}") if x),
                    time=f"{hd.get('TIMESYS')} {hd.get('TIMEREF')} BJDREFI={hd.get('BJDREFI')}", sector=hd.get("SECTOR"),
                    cols=",".join(h[1].data.columns.names[:14]))
    except Exception as e:  # noqa
        return dict(release="", time="", sector="", cols=repr(e))


USE = {
    "tess2020324010417-s0032-0000000004206066-0200-s_lc.fits": ("S32 SPOC 2-min light curve", "PDCSAP_FLUX (primary), SAP_FLUX (variant)", "lightkurve DEFAULT_BITMASK 17087", "SPOC optimal aperture (18 px); CROWDSAP 0.997, FLFRCSAP 0.922"),
    "tess2020324010417-s0032-0000000004206066-0200-s_tp.fits": ("S32 SPOC 2-min target pixel file", "FLUX (calibrated, background-subtracted)", "DEFAULT_BITMASK", "own circular apertures r=1-3 px; PRF fits"),
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_lc.fits": ("S6 TESS-SPOC 30-min FFI light curve", "PDCSAP_FLUX / SAP_FLUX", "DEFAULT_BITMASK", "TESS-SPOC optimal aperture 18 px; CROWDSAP 0.995"),
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0006_tess_v1_tp.fits": ("S6 TESS-SPOC target pixels", "FLUX", "DEFAULT_BITMASK", "aperture bit 2 used for TRICERATOPS/census"),
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_lc.fits": ("S32 TESS-SPOC 10-min light curve (comparison)", "PDCSAP_FLUX", "DEFAULT_BITMASK", "TESS-SPOC optimal aperture; CROWDSAP 0.993"),
    "hlsp_tess-spoc_tess_phot_0000000004206066-s0032_tess_v1_tp.fits": ("S32 TESS-SPOC target pixels (comparison)", "FLUX", "DEFAULT_BITMASK", "own apertures"),
    "hlsp_qlp_tess_ffi_s0006-0000000004206066_tess_v01_llc.fits": ("S6 QLP 30-min (comparison)", "KSPSAP_FLUX / SAP_FLUX", "DEFAULT_BITMASK + QLP bit 13 (4096)", "QLP best aperture"),
    "hlsp_qlp_tess_ffi_s0032-0000000004206066_tess_v01_llc.fits": ("S32 QLP 10-min (comparison)", "KSPSAP_FLUX / SAP_FLUX", "DEFAULT_BITMASK + bit 13", "QLP best aperture"),
    "hlsp_qlp_tess_ffi_s0098-0000000004206066_tess_v01_llc.fits": ("S98 QLP 200-s (TGLC-era, TICA v1.4.0)", "SYS_RM_FLUX (primary), SAP/DET (variants)", "frozen mask 1610619071 (bits 0-5,7,11,12,29,30) primary; variants: DEFAULT_BITMASK+2^30, +straylight", "QLP BESTAP '3x3 square' (TGLC PSF-based)"),
    "hlsp_gsfc-eleanor-lite_tess_ffi_s0006-0000000004206066_tess_v1.0_lc.fits": ("S6 eleanor-lite (comparison)", "CORR_FLUX", "QUALITY==0", "eleanor aperture"),
    "hlsp_tglc_tess_ffi_gaiaid-3220388198192519424-s0006-cam1-ccd1_tess_v1_llc.fits": ("S6 TGLC (comparison)", "cal_aper_flux (and cal_psf_flux variant)", "TESS_flags DEFAULT_BITMASK, TGLC_flags==0", "TGLC 3x3 aperture on ePSF-decontaminated image"),
    "hlsp_t16_tess_ffi_s0006-cam1-ccd1-03220388198192519424_tess_v01_llc.fits": ("S6 T16 (inventory only)", "-", "-", "-"),
    "hlsp_tasoc_tess_ffi_tic00004206066-s0006-cam1-ccd1-c1800_tess_v05_cbv-lc.fits": ("S6 TASOC (inventory only)", "-", "-", "-"),
    "hlsp_tasoc_tess_ffi_tic00004206066-s0006-cam1-ccd1-c1800_tess_v05_ens-lc.fits": ("S6 TASOC (inventory only)", "-", "-", "-"),
    "hlsp_tars_tess_ffi_s0006-0000000004206066_tess_v01_lc.fits": ("S6 TARS (timing check only: TIME not target-barycentred)", "-", "-", "-"),
    "hlsp_tars_tess_ffi_s0032-0000000004206066_tess_v01_lc.fits": ("S32 TARS (timing check only)", "-", "-", "-"),
}
for fn, (use, flux, qual, ap) in USE.items():
    p = os.path.join(AUD, "data", "mast", fn)
    fm = fits_meta(p)
    r = dl.get(fn, {})
    pv = last.get(os.path.join("data", "mast", fn), {})
    rows.append(dict(dataset=use, file=f"data/mast/{fn}", source=r.get("dataURI", ""), url=pv.get("url", r.get("url", "")),
                     release_version=fm["release"], download_utc=pv.get("utc_end", "cached (same session)"), sha256=sha256_file(p),
                     same_as_deep_followup=r.get("same_as_old"), units="relative flux after median normalisation (e-/s in file)",
                     time_system=fm["time"] + " (BTJD = BJD_TDB - 2457000)", flux_definition=flux, quality_filtering=qual, aperture=ap))
for f, c in tc.items():
    p = os.path.join(AUD, f)
    pv = [x for x in prov if x.get("dest", "").endswith(f"tesscut_s{c['sector']}.zip")]
    rows.append(dict(dataset=f"S{c['sector']} TESScut 21x21 FFI cutout (SPOC FFIs, astrocut {c['procver']})", file=f, source="MAST TESScut",
                     url=pv[-1]["url"] if pv else "", release_version=f"FFI_TYPE={c['ffi_type']}; cutout DATE={c['date']}",
                     download_utc=pv[-1]["utc_end"] if pv else "", sha256=sha256_file(p),
                     same_as_deep_followup=f"FLUX identical={c.get('flux_identical')}, TIME identical={c.get('time_identical')}",
                     units="e-/s per pixel (not background-subtracted)", time_system=f"{c['timesys']} BJDREFI={c['bjdrefi']}; barycentred to CCD centre (astrocut)",
                     flux_definition="own aperture photometry (r = 1-3 px) after per-cadence background (median of faint pixels > 5 px)",
                     quality_filtering="DEFAULT_BITMASK; straylight (2048) dropped for S98 in localisation", aperture="circles r=1.0-3.0 px; PRF fits within 4 px"))
for dest, p in sorted(last.items()):
    if dest.startswith("data/mast/") or "tesscut_s" in dest:
        continue
    use = ("Gaia DR3/DR2 (ESA TAP)" if dest.startswith("data/gaia") else "VizieR / MAST catalogue" if dest.startswith("data/catalogs")
           else "reference document / method source" if dest.startswith("data/refs") else "TESS PRF model (MAST start_s0004)" if dest.startswith("data/prf")
           else "novelty refresh" if dest.startswith("data/novelty") else "other")
    rows.append(dict(dataset=use, file=dest, source="", url=p.get("url", ""), release_version="", download_utc=p.get("utc_end", ""),
                     sha256=p.get("sha256", ""), same_as_deep_followup="", units="", time_system="", flux_definition="",
                     quality_filtering="", aperture=p.get("note", "")))
# re-used earlier-session files (not re-downloaded; hashes recorded here)
for f, use in [(os.path.join(OLD_DEEP, "queries", "asassn_skypatrol_v2_lc.csv"), "ASAS-SN Sky Patrol v2 light curve (re-used, pyasassn query 2026-09-25T03:20:56Z by deep follow-up)"),
               (os.path.join(OLD_COMP, "data", "way2026_overluminous_row.csv"), "Way et al. 2026 overluminosity row (re-used; Zenodo 18500082)"),
               (os.path.join(OLD_COMP, "data", "way2026_systems_row.csv"), "Way et al. 2026 systems row (re-used)"),
               (os.path.join(OLD_COMP, "sed_out", "best_fit_average.dat"), "ARIADNE SED result (re-used, not re-run)")]:
    rows.append(dict(dataset=use, file=os.path.relpath(f, AUD), source="earlier session", url="", release_version="", download_utc="",
                     sha256=sha256_file(f), same_as_deep_followup="re-used", units="", time_system="", flux_definition="", quality_filtering="", aperture=""))
keys = ["dataset", "file", "source", "url", "release_version", "download_utc", "sha256", "same_as_deep_followup", "units", "time_system",
        "flux_definition", "quality_filtering", "aperture"]
with open(os.path.join(AUD, "DATA_PROVENANCE.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "rows written")
