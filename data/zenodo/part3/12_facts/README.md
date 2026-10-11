# 12_facts: the authoritative fact table and the Note-to-package map

| file | content |
|---|---|
| `FINAL_FACTS.csv` | The reconciled fact table: 197 facts (192 FINAL, 2 UNRESOLVED, 3 NOT_AVAILABLE; includes `A.archival_FAP` added after the pass-3 review). Copied **unchanged** from the author's reconciliation folder. Columns: `fact_id`, `group`, `quantity`, `value`, `uncertainty`, `units`, `status`, printed form (`manuscript_form`), `evidence`, `code`, `rerun`, `supersedes`, `verification`, `topic`, `notes`. |
| `MANUSCRIPT_NUMBER_AUDIT.csv` | Every number or quantity printed in the Note (135 rows; pass-3 revision), with its fact id, final value and match type. Copied unchanged. |
| `NOTE_NUMBER_MAP.csv` | The same 135 rows plus the package file and key that hold each number (first mapped evidence file; others in `other_package_files`), and `check_type`: the audit match type, or `BOUND` when the printed value is an upper limit ("≲1%", "K<1.2") that the check accepts if the file holds the bounded value itself. The table in `../README.md` is built from this file. |
| `PATH_MAP.csv` | Original folder prefix → package folder prefix, for every copied source folder. |
| `EXCLUDED_FILES.csv` | Every file or folder in the copied source trees that is **not** in the package, with size, sha256 and the reason: raw FITS, MCMC chains, pickles, third-party papers and documentation, large tables, internal search bookkeeping. |
| `FACTS_IN_THIS_SECTION.csv` | The facts not assigned to a topic folder (manuscript metadata). |
| `COVERAGE_CHECK.json` | For every row of the Note-number table: the printed numbers and the package file in which a matching number was found. Method in `summary.method`. Result: 0 gaps. |

**Paths.** The `evidence`, `code` and `rerun` columns of `FINAL_FACTS.csv` give paths relative to the author's project tree (`addenda/...`, `candidates/...`, `config/...`). Some are absolute paths on the author's machine. They are kept as written. To find a file in this package, use `PATH_MAP.csv`, or the `evidence_in_package` column of each folder's `FACTS_IN_THIS_SECTION.csv`, which already translates the paths.

Paths under `candidates/`, `config/`, `data/` and `logs/` belong to the author's wider TESS search. They are not included.
