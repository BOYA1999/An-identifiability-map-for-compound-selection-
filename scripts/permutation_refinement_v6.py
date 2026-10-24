"""Refinement of the scaffold-missingness permutation statistic.

Adds a higher-precision permutation null (20,000 draws) with an explicitly defined
statistic, null hypothesis and one-sided empirical p-value, without altering any
previously verified artifact.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'artifacts/expansionrx_v6_strengthening_v1'
STRUCT = ROOT / 'artifacts/expansionrx_structure_audit/structure_membership.csv'
RAW = ROOT / 'data/expansionrx_openadmet/expansion_data_raw.csv'
SCORES = ROOT / 'artifacts/expansionrx_later_fixed_scores_v2/policy_scores.csv'
N_DRAWS = 20000
SEED = 42

scores = pd.read_csv(SCORES)
names = scores['Molecule Name'].to_numpy()
raw = pd.read_csv(RAW, dtype=str, keep_default_na=False).set_index('Molecule Name').loc[names].reset_index()
struct = pd.read_csv(STRUCT)
struct = struct[struct.official_split == 'test'].set_index('Molecule Name').loc[names].reset_index()
assert struct['Molecule Name'].tolist() == list(names)

missing = (raw['HLM CLint'] == '').to_numpy()
group = struct['murcko_nonstereo'].fillna('').to_numpy()
keep = group != ''
codes, uniques = pd.factorize(group[keep])
miss = missing[keep]
n_groups = len(uniques)
counts = np.bincount(codes, minlength=n_groups)
group_has5 = counts >= 5

# Statistic T: share of all missing clearance records that lie in scaffold groups of size >= 5
# whose within-group missing rate is at least 0.8.
observed_missing = int(miss.sum())


def statistic(perm_missing):
    sums = np.bincount(codes, weights=perm_missing, minlength=n_groups)
    with np.errstate(divide='ignore', invalid='ignore'):
        rates = np.where(counts > 0, sums / counts, 0.0)
    ok_groups = group_has5 & (rates >= 0.8)
    return float(perm_missing[ok_groups[codes]].sum() / observed_missing)

T_obs = statistic(miss.astype(float))
rng = np.random.RandomState(SEED)
draws = np.empty(N_DRAWS, dtype=float)
base = miss.astype(int)
for i in range(N_DRAWS):
    draws[i] = statistic(rng.permutation(base).astype(float))
p_one_sided = (1.0 + int((draws >= T_obs).sum())) / (1.0 + N_DRAWS)
result = {
    'statistic': 'share of missing clearance records in Murcko scaffold groups with >=5 records and within-group missing rate >= 0.8',
    'null_hypothesis': 'clearance missingness is exchangeable across the 2,282 later compounds given the observed scaffold assignment',
    'alternative': 'missingness is more concentrated in high-missing scaffold groups than exchangeability predicts (one-sided)',
    'draws': N_DRAWS, 'seed': SEED, 'rng': 'numpy.random.RandomState',
    'n_compounds_with_scaffold': int(keep.sum()), 'n_groups': int(n_groups),
    'groups_size_at_least_5': int(group_has5.sum()), 'n_missing_records': observed_missing,
    'observed_statistic': T_obs,
    'permutation_mean': float(draws.mean()), 'permutation_sd': float(draws.std(ddof=1)),
    'permutation_p95': float(np.percentile(draws, 95)), 'permutation_p99': float(np.percentile(draws, 99)),
    'permutation_p999': float(np.percentile(draws, 99.9)),
    'permutation_max': float(draws.max()),
    'one_sided_empirical_p': p_one_sided,
    'p_value_upper_bound': 1.0 / (1.0 + N_DRAWS),
    'inputs_sha256': {'structure_membership': hashlib.sha256(STRUCT.read_bytes()).hexdigest(),
                      'later_raw': hashlib.sha256(RAW.read_bytes()).hexdigest(),
                      'later_scores': hashlib.sha256(SCORES.read_bytes()).hexdigest()},
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(OUT / 'permutation_refinement.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
pd.DataFrame({'draw': np.arange(N_DRAWS), 'statistic': draws}).to_csv(OUT / 'permutation_refinement_draws.csv', index=False)
print(json.dumps(result, indent=2))
