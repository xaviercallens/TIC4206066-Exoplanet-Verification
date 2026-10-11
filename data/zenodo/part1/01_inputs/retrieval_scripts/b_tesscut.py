"""(b) TESScut FFI cutouts, sectors 6/32/98, 21x21 and 31x31, via astroquery.mast.Tesscut."""
import os, sys, json, zipfile, time
sys.path.insert(0, os.path.dirname(__file__))
from common import *
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.mast import Tesscut
import astroquery

c = SkyCoord(RA, DEC, unit="deg")
res = []
for size in (21, 31):
    for sector in (6, 32, 98):
        d = os.path.join(ROOT, "tesscut", f"{size}x{size}")
        os.makedirs(d, exist_ok=True)
        # the exact URL astroquery 0.4.11 builds (cutouts.py: REQUEST_URL + 'astrocut?' + params)
        url = f"{Tesscut._service_api_connection.REQUEST_URL}astrocut?y={size}&x={size}&sector={sector}&ra={c.ra.deg}&dec={c.dec.deg}"
        for attempt in range(4):
            t0 = utcnow()
            try:
                tab = Tesscut.download_cutouts(coordinates=c, size=size, sector=sector, path=d, inflate=False)
                zp = str(tab["Local Path"][0])
                break
            except Exception as ex:
                print("retry", size, sector, repr(ex), flush=True)
                time.sleep(10 * (attempt + 1))
        else:
            res.append(dict(size=size, sector=sector, error="download failed", url=url)); continue
        t1 = utcnow()
        zsha, zbytes = sha256_file(zp), os.path.getsize(zp)
        with zipfile.ZipFile(zp) as z:
            names = z.namelist()
            z.extractall(d, members=names)
        for n in names:
            fp = os.path.join(d, n)
            e = dict(relpath=os.path.relpath(fp, ROOT), url=url, utc_start=t0, utc=t1, bytes=os.path.getsize(fp),
                     sha256=sha256_file(fp), zip_sha256=zsha, zip_bytes=zbytes, tool=f"astroquery {astroquery.__version__} Tesscut.download_cutouts",
                     size=size, sector=sector)
            log_prov(e); res.append(e)
            print(json.dumps(e), flush=True)
        os.remove(zp)  # zip sha256 is recorded; keep only the inflated FITS to save shared disk
json.dump(res, open(os.path.join(ROOT, "tesscut", "tesscut_download_log.json"), "w"), indent=1)
print("DONE", flush=True)
