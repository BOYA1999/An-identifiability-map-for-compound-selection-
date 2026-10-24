# An identifiability map for compound selection

Analysis code for **What a released ADMET dataset can decide about compound selection: an identifiability map**.

This repository reproduces the endpoint-profile maps, threshold sensitivity, identification distances, conditional coverage requirements, feasible greedy measurement frontier, sharpness witnesses, missingness analyses, paired-bound decomposition and figures from the study's frozen policy records. It uses the same scores and selections; it does not retrain the prediction models or generate new predictions.

## Run

Use Python 3.12 and a new working directory:

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

The runner verifies input SHA-256 checksums, downloads the pinned original data, regenerates the analyses and figures, and checks 27 numerical outputs against the frozen reference. It writes results to `work/`, figures to `work/reporting/figures/`, and the verification result to `work/verification.json`. Existing working directories are not overwritten. To run again, choose a new directory with `--workdir another_run`.

The comparison requires exact CSV bytes; JSON scientific values use an absolute tolerance of 1e-12. Generation timestamps and producer-script identity fields are excluded from the JSON comparison. Graphical appearance can vary with installed fonts; figure images are not compared byte for byte.

## Inputs and sources

- Original dataset: [OpenADMET–ExpansionRx, pinned release](https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/tree/6b898ccc43d10d25b230fb09e22a6e30c30022b5).
- Frozen inputs: 13 CSV files containing the saved scores, selections, outcome records, membership mappings and baseline summaries required by the analysis. These author-generated records are not part of the original dataset. The runner automatically downloads the [frozen input archive](https://github.com/BOYA1999/An-identifiability-map-for-compound-selection-/releases/download/reproducibility-20261024/frozen_inputs.zip). To use a previously downloaded copy, run `python reproduce.py --inputs /path/to/frozen_inputs.zip`.
- Exact download addresses, archive contents and checksums are recorded in `data_sources.json`.

The original dataset README declares CC BY 4.0. The repository's existing [MIT license](LICENSE) covers code; it does not replace the original data provider's terms or assign a new licence to separately supplied records.

## Code and outputs

`scripts/` contains the profile, measurement and permutation analyses plus the observation-scope plotting code. `reporting/analysis/` contains the paired decomposition and remaining figure code. `reproduce.py` prepares the inputs and runs these scripts in a separate directory. `expected_results.json` contains only the hashes and reference values needed to check the regenerated outputs.

The figure order follows the current Digital Discovery manuscript: Figure 1, observation scope; Figure 2, identifiability map; Figure 3, internal positive control and conditional measurement analyses; Figure 4, evidence gaps. Supplementary Figures S1–S4 are also regenerated. The `v5` and `v6` names inside the code identify the frozen analysis stages, not alternative manuscript versions.

The known-success difference and policy-exclusive unknown outcomes jointly determine identification. The 91.9% figure concerns a fixed-difference, unknown-as-failure route. The 217-measurement schedule is feasible for the 44 nonzero-difference comparisons when all resolved outcomes fail; it is not a proven minimum or an expected acquisition cost. Reproducing these calculations is not external scientific validation.
