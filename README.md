# An identifiability map for compound selection

Code for **What a released ADMET dataset can decide about compound selection: an identifiability map**.

The analysis asks when incomplete ADMET measurements are sufficient to compare two sets of selected compounds. The scripts use the study's saved model scores and selections to reproduce the tables and figures.

## Run

With Python 3.12 installed, run these commands from the repository directory:

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

The script downloads the required data, checks the files, and runs the analyses. Results are saved in `work/`, figures in `work/reporting/figures/`, and the results of 27 original output checks in `work/verification.json`. The additional decision, composition and measurement replay results are saved in `work/analysis_depth/`, with eight comparisons against the recorded outputs in `work/analysis_depth/verification.json`.

Each run needs a new output directory. To run the analysis again, use `python reproduce.py --workdir another_run`.

## Data

- [OpenADMET–ExpansionRx](https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/tree/6b898ccc43d10d25b230fb09e22a6e30c30022b5): the dataset version used in the study.
- [Analysis inputs](https://github.com/BOYA1999/An-identifiability-map-for-compound-selection-/releases/tag/reproducibility-20261024-r1): 13 CSV files with the scores, selections, outcome records, group assignments and summary tables prepared for this study.
- `data_sources.json`: download addresses, file lists and checksums.

Both sources are downloaded automatically. To use an input archive you have already downloaded, run `python reproduce.py --inputs /path/to/inputs.zip`. If the source CSVs are already available, `python reproduce.py --inputs /path/to/inputs.zip --source-dir /path/to/source-files` uses the two local source files after checking their checksums.

## Files

`scripts/` contains the endpoint-profile, measurement and permutation analyses. `depth/` contains the joint policy feasibility, worst-case shortfall, endpoint-composition and internal measurement-replay analyses, together with the decision figure. `reporting/analysis/` contains the paired comparisons and plotting code. `reproduce.py` runs the workflow, and `expected_results.json` provides the reference outputs used to check the results.

The original CSV outputs are checked by file checksum. Numerical values in JSON outputs are compared with an absolute tolerance of 1e-12; timestamps and script labels are omitted from this comparison. Additional analysis outputs are checked by checksums of selected fields, including all 3,000 internal replay rows. Feasible-solution witnesses may vary across solver versions; the check compares feasibility rather than the particular completion. Figure appearance may vary with the fonts installed on your computer. The measurement replay uses previously observed development records, not newly measured outcomes or external validation.

## License

The code is available under the [MIT license](LICENSE). The original OpenADMET–ExpansionRx dataset README states CC BY 4.0; please refer to the source for its data-use terms.
