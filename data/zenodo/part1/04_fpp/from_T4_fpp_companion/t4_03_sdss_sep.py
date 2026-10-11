"""C020: separation of the extra SDSS PhotoObjAll entry from the saturated star (read-only raw SkyServer JSON)."""
import json, glob, os
from astropy.coordinates import SkyCoord
import astropy.units as u
A = "/Users/paulsmith/space <3/tess-jev/addenda/archival_extra_tic4206066"
res = {}
for f in sorted(glob.glob(A + "/spectra_activity/raw/sdss_dr1*_photo.json")):
    rows = json.load(open(f))[0]["Rows"]
    s = [SkyCoord(r["ra"] * u.deg, r["dec"] * u.deg) for r in rows]
    res[os.path.basename(f)] = dict(
        star_objid=rows[0]["objid"], extra_objid=rows[1]["objid"],
        sep_star_to_extra_arcsec=round(s[0].separation(s[1]).arcsec, 3),
        query_distance_arcsec=[round(r["distance"] * 60, 3) for r in rows],
        extra_psfMag_r_i=[rows[1].get("psfMag_r"), rows[1].get("psfMag_i")], keys=list(rows[1].keys()))
json.dump(res, open(os.path.dirname(os.path.abspath(__file__)) + "/t4_03_sdss_sep.json", "w"), indent=1)
print(json.dumps(res, indent=1))
