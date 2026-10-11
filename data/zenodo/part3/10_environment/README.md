# 10_environment: software used

| file | content |
|---|---|
| `VERSIONS.txt` | OS (macOS 26.3.1, arm64), Python versions, key package versions per environment, and the LaTeX/counting tools: Tectonic 0.17.0, AASTeX 7.0.2, TeXcount 3.1.1, pdftotext 26.08.0. |
| `pip_freeze_main_venv.txt` | `python -m pip freeze` of the main project environment (Python 3.12.13). |
| `pip_freeze_trice_env.txt` | Package list of the TRICERATOPS environment (Python 3.11.15; triceratops 1.1.0, numpy 1.26.4). pip is not installed in this env, so the list comes from `importlib.metadata`; the format is the same as pip freeze. |
| `pip_freeze_audit_env.txt` | Package list of the opus55 audit environment (Python 3.12.13), taken the same way. |

Other environment records inside copied folders:
- `../13_supporting_pipelines/PUBLICATION_AUDIT/logs/audit_env_freeze.txt`;
- `../06_pixels_localisation/from_tic4206066_pixel_recheck/env/pip_freeze.txt`.

The astroARIADNE 1.5.0 environment used for the SED fit was not recorded in full (see `VERSIONS.txt`).
