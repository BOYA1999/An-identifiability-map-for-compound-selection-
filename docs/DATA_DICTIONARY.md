# Data dictionary

The repository is the verbatim merged payload of the three ESI archives (AF2 data/code,
AF3 and AF4 fitted assets). Paths below are relative to `package/ExpansionRx_reproducibility/`.

## Top-level directories

| Path | Contents |
|---|---|
| `data/` | Released OpenADMET–ExpansionRx source tables and README; curation ledger |
| `artifacts/` | Saved selections, sealed scores, fitted model assets, pilot failures and audit artifacts |
| `protocol/` | Frozen scientific protocols (scenario, policies, thresholds, evaluation) |
| `reporting/` | Derived reporting tables, figures and source manifests used in the paper |
| `scripts/` | Producer, replay and reproducibility scripts |
| `verification/` | Independent-check records (23/23, 8/8, 18/18) and their scripts |
| `provenance/` | Historical frozen copies, including the original v6 refinement record |
| `licenses/` | Bundled third-party license (Chemprop, MIT) |
| `environment.json` | Recorded Python/scientific stack (training stack; not needed for replay) |
| `MATERIALS_MANIFEST.json` | Exact combined-file manifest: path, bytes, SHA-256 |
| `README.md` | Original package README with the five-command verified route |

## Reporting data index (paper ESI data S1–S33)

| Data | File in `reporting/tables/` (unless noted) | Data rows |
|---|---|---|
| S1 | `supplementary_data_s1_development_policy_metrics.csv` | 30 |
| S2 | `supplementary_data_s2_development_pairwise_identification.csv` | 135 |
| S3 | `supplementary_data_s3_development_selection_membership.csv` | 4,625 |
| S4 | `supplementary_data_s4_later_policy_metrics.csv` | 30 |
| S5 | `supplementary_data_s5_later_pairwise_identification.csv` | 135 |
| S6 | `supplementary_data_s6_later_selection_membership.csv` | 8,190 |
| S7 | `supplementary_data_s7_later_probability_promises.csv` | 18 |
| S8 | `supplementary_data_s8_later_scored_outcomes.csv` | 2,282 |
| S9 | `supplementary_data_s9_later_policy_scores.csv` | 2,282 |
| S10 | `supplementary_data_s10_later_rankings_and_selection.csv` | 22,820 |
| S11 | `supplementary_data_s11_structure_membership.csv` | 7,608 |
| S12 | `supplementary_data_s12_calibration_development_scores.csv` | 2,390 |
| S13 | `paired_decomposition_all_135.csv` | 135 |
| S14 | `joint_unknown_endpoint_patterns.csv` | 6 |
| S15 | `representative_pairs.csv` | 3 |
| S16 | `threshold_source_mapping.csv` | 3 |
| S17 | `supplementary_data_s17_development_prediction_metrics.csv` | 6 |
| S18 | `profile_identifiability.csv` | 21 |
| S19 | `profile_outcome_counts.csv` | 7 |
| S20 | `threshold_sensitivity.csv` | 10 |
| S21 | `identification_distance_primary45.csv` | 45 |
| S22 | `selected_unknown_rates.csv` | 30 |
| S23 | `profile_pair_intervals.csv` | 945 |
| S24 | `threshold_pair_intervals.csv` | 1,350 |
| S25 | `profile_map_later.csv` | 7 |
| S26 | `profile_map_development.csv` | 7 |
| S27 | `coverage_requirement_primary45.csv` | 45 |
| S28 | `measurement_frontier.csv` | 217 |
| S29 | `sharpness_witnesses.csv` | 270 |
| S30 | `missingness_structure_summary.csv` | 1 |
| S31 | `naive_report_magnitudes.csv` | 6 |
| S32 | `permutation_refinement.json` | — |
| S33 | `permutation_refinement_draws.csv` | 20,000 |
