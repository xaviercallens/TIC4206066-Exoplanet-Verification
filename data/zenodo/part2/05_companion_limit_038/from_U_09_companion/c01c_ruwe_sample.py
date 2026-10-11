"""Audit 9: RUWE / AEN / RV comparison sample from the Gaia archive using the indexed random_index column."""
import hashlib, json, datetime
from pathlib import Path
from astroquery.gaia import Gaia
HERE = Path(__file__).resolve().parent; D = HERE / "data"
q = ("SELECT source_id, phot_g_mean_mag, bp_rp, parallax, ruwe, astrometric_excess_noise, astrometric_excess_noise_sig, "
     "ipd_gof_harmonic_amplitude, ipd_frac_multi_peak, non_single_star, radial_velocity_error, rv_nb_transits, "
     "rv_amplitude_robust, rv_chisq_pvalue, rv_renormalised_gof, grvs_mag, rv_template_teff FROM gaiadr3.gaia_source "
     "WHERE random_index < 300000000 AND phot_g_mean_mag BETWEEN 10.0 AND 11.0 AND bp_rp BETWEEN 1.55 AND 1.95 AND parallax > 5")
Gaia.ROW_LIMIT = -1
t = Gaia.launch_job_async(q).get_results()
p = D / "ruwe_comparison_random.ecsv"; t.write(p, format="ascii.ecsv", overwrite=True)
open(D / "PROVENANCE.jsonl", "a").write(json.dumps(dict(name="ruwe_comparison_random", url="https://gea.esac.esa.int/tap-server/tap", query=q, file=p.name, nrows=len(t), sha256=hashlib.sha256(p.read_bytes()).hexdigest(), utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))) + "\n")
print(len(t))
