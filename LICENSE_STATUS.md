# License and redistribution status

This repository is prepared for public release after the manuscript is accepted, but the final licenses
have **not** been chosen yet. Do not assume an open-source license for the added code or fitted models.

## Components

| Component | Status |
|---|---|
| OpenADMET–ExpansionRx source data (`package/.../data/`) | The archived release README declares CC BY 4.0. The separately linked LICENSE.md was unavailable at the source audit, so redistribution of these files should be confirmed against the release's current terms before publishing this repository. |
| Chemprop package license | Included verbatim: `package/ExpansionRx_reproducibility/licenses/chemprop_LICENSE.txt` (MIT, Chemprop Development Team). |
| Newly generated analysis code (`scripts/`, `reporting/analysis/`, `verification/`) | No final license yet. The authors should add a code license (e.g., MIT or Apache-2.0) and a copyright line. |
| Fitted model assets (`artifacts/**/*.joblib`, checkpoints) | No final license yet. The authors should confirm they have the right to redistribute them and add a license. |
| Derived tables and figures (`reporting/tables/`, figures) | No final license yet; the authors should choose, e.g., CC BY 4.0, consistent with the source data. |

## Before making the repository public

1. Choose and add `LICENSE` (code), `LICENSE-DATA` (data/tables) and, if needed, a model-weights note.
2. Confirm redistribution rights for the source data and the nine fitted assets.
3. Add the final author list, paper DOI and repository DOI to README and CITATION.
4. Optionally archive a release on Zenodo to obtain a DOI.

Until these steps are complete, the repository should remain private.
