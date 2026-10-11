"""The Lean certification step must fail closed: no Lean run => not certified."""
import os
import shutil

import pytest

from scripts.run_tic4206066_workflow import TIC4206066Workflow


def test_lean_gate_fails_closed_without_lake(monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    res = TIC4206066Workflow().step8_lean4_formal_verification()
    assert res.passed is False
    assert res.details["lean_elaboration_clean"] is False
    assert res.details["lean_elaboration_exit_code"] == -1


def test_lean_source_has_no_sorry_and_six_theorems():
    src = open(os.path.join("formal", "ANSE", "ExoplanetTransit.lean"), encoding="utf-8").read()
    code = "\n".join(l.split("--")[0] for l in src.splitlines())
    assert "sorry" not in code
    assert code.count("\ntheorem ") == 6
