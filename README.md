# What a released ADMET dataset can decide — reproducibility package

Code, saved results and fitted assets for the manuscript **"What a released ADMET dataset can decide about compound selection: an identifiability map"** (submitted to *Digital Discovery*).

The study asks whether a public, partially observed ADMET release can identify differences between fixed compound-selection policies. Using the OpenADMET–ExpansionRx later cohort (2,282 compounds, ten fixed policies, three budgets, 45 pairwise comparisons) it computes sharp finite-pool identification bounds and maps what the released records decide, where the deciding power sits, and which records would have to be measured to change an answer.

**This repository is a reproducibility package, not a re-training release.** All analyses reuse frozen scores and selections; no model is retrained and no prediction is regenerated. The verified route below reproduces the saved results and the reported derived tables in minutes.

## What is included

| Part | Contents |
|---|---|
| `package/ExpansionRx_reproducibility/data/` | Released source tables and curation record |
| `package/ExpansionRx_reproducibility/reporting/` | Derived reporting tables (ESI data S1–S33), figure sources and manifests |
| `package/ExpansionRx_reproducibility/scripts/` | Analysis producers, portable saved-result replay, package verifier |
| `package/ExpansionRx_reproducibility/verification/` | Independent-check records: 23/23, 8/8 and 18/18 checks |
| `package/ExpansionRx_reproducibility/artifacts/` | Saved selections, sealed scores and the nine fitted model assets (from ESI AF3–AF4) |
| `checksums/SHA256SUMS` | SHA-256 of every packaged file (matches the internal manifest) |

## Quick start

Requirements: Python 3.12. The saved-result replay uses the Python standard library only.

```bash
cd package/ExpansionRx_reproducibility
python -B -X utf8 verify_package.py     # exact file set, sizes and SHA-256 (335 files)
python -B -X utf8 replay_saved.py       # re-derives saved selections, labels and 135 contrasts
```

`verify_package.py` enforces the complete combined payload (AF2 + AF3 + AF4). It accepts only the
self-referential root manifest and the declared runtime directories as exceptions.

To rebuild the derived reporting tables (needs NumPy and pandas):

```bash
pip install -r requirements-derived.txt
python -B -X utf8 reproduce_derived.py  # fresh work copy, runs v5/v6 producers and checkers
```

See [docs/REPRODUCING.md](docs/REPRODUCING.md) for expected outputs, the three-checker record and troubleshooting.

## Verification status

| Check | Result |
|---|---|
| Combined payload manifest | 335 files, exact sizes and SHA-256 |
| Independent identification check (v5) | 23/23 passed |
| Independent strengthening check (v6) | 8/8 passed |
| Independent permutation-refinement check | 18/18 passed |
| Portable replay | saved selections, labels and 135 contrasts reproduced |

These are implementation-independence checks on the same frozen inputs. They are **not** external
scientific validation, new experiments or an independent re-measurement of any assay.

## Data provenance and licensing — please read

- The OpenADMET–ExpansionRx source data come from the public release pinned at revision
  `6b898ccc43d10d25b230fb09e22a6e30c30022b5`; the archived release README declares CC BY 4.0.
- The bundled Chemprop license is included at `package/ExpansionRx_reproducibility/licenses/chemprop_LICENSE.txt` (MIT).
- **Newly generated code, fitted models and derived tables have no final license yet.** The authors must
  choose and add code and data licenses, and confirm redistribution rights, before this repository is made
  public. See [LICENSE_STATUS.md](LICENSE_STATUS.md). Until then, treat the content as "all rights reserved
  for review purposes".

## Citing

The manuscript is under submission; a DOI and author list will be added on acceptance. A citation
template is provided in [CITATION.cff.template](CITATION.cff.template). Until publication, please cite the
manuscript title and the pinned data revision above.

## Repository status

- Author names, affiliations and the public repository URL/DOI are pending and are marked in the
  submission materials.
- `package/` is the verified combined reproducibility payload (335 files, each bound by
  `MATERIALS_MANIFEST.json`); it is unchanged from the manuscript's submission payload. The wrapper
  (README, docs, CI, checksums, license notes, Makefile) is added only for GitHub hosting.
- `package/ExpansionRx_reproducibility/README.md` is the historical package README and is byte-bound by
  the manifest (and therefore intentionally not edited). The current manuscript, scope and commands are
  described in this top-level README and in `docs/REPRODUCING.md`.
