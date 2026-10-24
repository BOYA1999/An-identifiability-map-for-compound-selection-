# Reproducing the results

## 1. Environment

- Python 3.12.
- `verify_package.py` and `replay_saved.py`: standard library only.
- `reproduce_derived.py`: NumPy and pandas (see `requirements-derived.txt`).
- Plot regeneration: Matplotlib and Pillow.
- The original training stack is recorded in `environment.json`; **no training is needed or claimed** by
  any command here.

## 2. Verify the package

```bash
cd package/ExpansionRx_reproducibility
python -B -X utf8 verify_package.py
```

Expected: `{"status": "PASS", "immutable_files": 335, ...}`. The verifier checks the exact combined file
set, byte sizes and SHA-256 from `MATERIALS_MANIFEST.json`; unknown additions are rejected.

## 3. Replay the saved results

```bash
python -B -X utf8 replay_saved.py
```

The replay reads the released endpoint strings and sealed scores, refuses to run into an existing
`replay_output/`, and recomputes labels, saved selections and the 135 pairwise contrasts without fitting,
inference or model loading. Its scientific outputs are compared with the original reference bytes;
`input_bindings.json` records the curation and therefore has a different identity from the original
private closure.

## 4. Rebuild the derived analyses

```bash
python -B -X utf8 reproduce_derived.py
```

This refuses to run into an existing `regeneration_output/`, makes a separate work copy, moves the
packaged v5/v6 outputs into `regeneration_output/authoritative/`, then executes the v5 producer, the v6
producer and the original seeded 20,000-draw permutation refinement. It runs the three independent
checkers in the work copy and writes `regeneration_output/comparison.json`. The extracted immutable
package is preserved.

Expected check families: 23/23 (v5 identification), 8/8 (v6 strengthening), 18/18 (permutation refinement).

## 5. Figures

Figure sources for the manuscript are in `reporting/figures/`; plotting scripts are in
`reporting/analysis/` and need Matplotlib/Pillow.

## Troubleshooting

- Run commands from `package/ExpansionRx_reproducibility/`; paths are relative to that root.
- If a run complains about an existing `replay_output/` or `regeneration_output/`, move or delete the
  directory and rerun.
- On Windows, use `python -B -X utf8` to avoid encoding differences in CSV round-trips.
- `MATERIALS_MANIFEST.json` is self-referential and excludes itself plus the declared runtime directories.
