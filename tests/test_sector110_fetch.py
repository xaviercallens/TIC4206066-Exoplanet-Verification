"""Sector 110 fetch/evaluate must be blind, deterministic and fail closed."""
import json

import scripts.fetch_sector110 as f


def _p(name):
    return {"filename": name, "uri": "mast:x/" + name, "size": 1, "collection": "HLSP"}


def test_pipeline_choice_is_frozen_not_alphabetical():
    prods = [_p("hlsp_tars_x_lc.fits"), _p("hlsp_qlp_x_llc.fits"), _p("hlsp_tess-spoc_x_lc.fits"),
             _p("hlsp_eleanor-lite_x_lc.fits")]
    assert f.select_product(prods)["filename"] == "hlsp_tess-spoc_x_lc.fits"
    assert f.select_product(prods[:2])["filename"] == "hlsp_qlp_x_llc.fits"  # falls back by priority
    assert f.select_product([_p("hlsp_tars_x_lc.fits"), _p("hlsp_eleanor-lite_x_lc.fits")]) is None


def test_fetch_awaiting_when_nothing_published(monkeypatch):
    monkeypatch.setattr(f, "list_products", lambda sector: [])
    assert f.fetch(110) == f.EXIT_AWAITING


def test_evaluate_requires_manifest(monkeypatch, tmp_path):
    monkeypatch.setattr(f, "MANIFEST", str(tmp_path / "none.json"))
    assert f.evaluate() == f.EXIT_AWAITING


def test_evaluate_detects_tampered_input(monkeypatch, tmp_path):
    fits = tmp_path / "s0110.fits"
    fits.write_bytes(b"original")
    man = tmp_path / "m.json"
    man.write_text(json.dumps({"files": [{"path": str(fits), "sha256": f.sha256_file(str(fits))}]}))
    fits.write_bytes(b"tampered")
    monkeypatch.setattr(f, "MANIFEST", str(man))
    assert f.evaluate() == 3
