"""Sector 110 fetch/evaluate must be blind, deterministic and fail closed."""
import json

import scripts.fetch_sector110 as f


def _p(name):
    return {"filename": name, "uri": "mast:x/" + name, "size": 1, "collection": "HLSP"}


def test_only_spoc_2min_primary_is_selectable():
    spoc = "tess2026365000000-s0110-0000000004206066-0300-s_lc.fits"
    prods = [_p("hlsp_tars_x_lc.fits"), _p("hlsp_qlp_x_llc.fits"), _p("hlsp_tess-spoc_x_lc.fits"), _p(spoc)]
    assert f.select_product(prods)["filename"] == spoc
    # QLP / TESS-SPOC fallbacks apply only after 2027-02-28 per the pre-registration: never selectable here
    assert f.select_product(prods[:3]) is None


def test_fetch_awaiting_when_nothing_published(monkeypatch):
    monkeypatch.setattr(f, "list_products", lambda sector: [])
    assert f.fetch(110) == f.EXIT_AWAITING


def test_evaluate_refuses_unregistered_rules():
    assert f.evaluate() == 4
