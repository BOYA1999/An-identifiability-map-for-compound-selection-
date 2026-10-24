# An identifiability map for compound selection

Code for **What a released ADMET dataset can decide about compound selection: an identifiability map**.

The analysis asks when incomplete ADMET measurements are sufficient to compare two sets of selected compounds. The scripts use the study's saved model scores and selections to reproduce the tables and figures.

## Run

With Python 3.12 installed, run these commands from the repository directory:

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

The script downloads the required data, checks the files, and runs the analyses. Results are saved in `work/`, figures in `work/reporting/figures/`, and the results of 27 output checks in `work/verification.json`.

Each run needs a new output directory. To run the analysis again, use `python reproduce.py --workdir another_run`.

## Data

- [OpenADMET–ExpansionRx](https://huggingface.co/datasets/openadmet/openadmet-expansionrx-challenge-data/tree/6b898ccc43d10d25b230fb09e22a6e30c30022b5): the dataset version used in the study.
- [Analysis inputs](https://github.com/BOYA1999/An-identifiability-map-for-compound-selection-/releases/tag/reproducibility-20261024-r1): 13 CSV files with the scores, selections, outcome records, group assignments and summary tables prepared for this study.
- `data_sources.json`: download addresses, file lists and checksums.

Both sources are downloaded automatically. To use an input archive you have already downloaded, run `python reproduce.py --inputs /path/to/inputs.zip`.

## Files

`scripts/` contains the endpoint-profile, measurement and permutation analyses. `reporting/analysis/` contains the paired comparisons and plotting code. `reproduce.py` runs the workflow, and `expected_results.json` provides the reference outputs used to check the results.

CSV outputs are checked by file checksum. Numerical values in JSON outputs are compared with an absolute tolerance of 1e-12; timestamps and script labels are omitted from this comparison. Figure appearance may vary with the fonts installed on your computer.

## License

The code is available under the [MIT license](LICENSE). The original OpenADMET–ExpansionRx dataset README states CC BY 4.0; please refer to the source for its data-use terms.
