"""Development-cohort comparisons and measurement analyses.

Uses the study records and reported endpoint measurements to calculate:
  A. Endpoint-profile identifiability map in the independent internal development cohort
     (a positive control: does the same procedure ever identify comparisons?).
  B. Coverage requirement per comparison (how complete the opposite-side exclusive set must be).
  C. Measurement-value frontier: greedy best-case plan and the number of comparisons it can make
     identifiable for m resolved exclusive unknowns.
  D. Explicit sharpness witnesses: a completion attaining each bound endpoint, for all pairs/budgets.
  E. Structural clustering of clearance missingness (within-scaffold concentration, permutation null)
     and whether selected batches preferentially come from high-missingness scaffold groups.
  F. Apparent effect sizes of the two complete-case reports (how misleading a naive reading is).
"""
import hashlib
import itertools
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/expansionrx_openadmet'
SCORES = ROOT / 'artifacts/expansionrx_later_fixed_scores_v2'
V4OUT = ROOT / 'artifacts/expansionrx_later_identification_v1'
DEV = ROOT / 'artifacts/expansionrx_selection_identification_development_v1'
TRAIN = ROOT / 'artifacts/expansionrx_decision_support_v1/adjudicated_training.csv'
STRUCT = ROOT / 'artifacts/expansionrx_structure_audit/structure_membership.csv'
OUT = ROOT / 'artifacts/expansionrx_v6_strengthening_v1'
OUT.mkdir(exist_ok=False)

ENDPOINTS = ['LogD', 'KSOL', 'HLM CLint']
BUDGETS = [0.05, 0.10, 0.20]
PRIMARY = dict(logd_lo=1.0, logd_hi=3.0, ksol_min=50.0, hlm_max=10.0)
NUM_RE = re.compile(r'\s*(<=|>=|<|>|=)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*')
SHA = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def property_state(value, endpoint, logd_lo, logd_hi, ksol_min, hlm_max):
    if value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip():
        return -1
    m = NUM_RE.fullmatch(str(value))
    if m is None:
        raise ValueError(f'unparseable endpoint string: {value!r} for {endpoint}')
    op, num = m.groups()
    v = float(num)
    if endpoint != 'LogD' and (v < 0 or (v == 0 and op in ('<', '<='))):
        return -1
    if op in (None, '='):
        if endpoint == 'LogD':
            return int(logd_lo <= v <= logd_hi)
        if endpoint == 'KSOL':
            return int(v > ksol_min)
        return int(v < hlm_max)
    if endpoint == 'LogD':
        if op == '>':
            return 0 if v >= logd_hi else -1
        if op == '>=':
            return 0 if v > logd_hi else -1
        if op == '<':
            return 0 if v <= logd_lo else -1
        if op == '<=':
            return 0 if v < logd_lo else -1
    if endpoint == 'KSOL':
        if op == '>':
            return 1 if v >= ksol_min else -1
        if op == '>=':
            return 1 if v > ksol_min else -1
        if op in ('<', '<='):
            return 0 if v <= ksol_min else -1
    if op == '<':
        return 1 if v <= hlm_max else -1
    if op == '<=':
        return 1 if v < hlm_max else -1
    if op in ('>', '>='):
        return 0 if v >= hlm_max else -1
    raise AssertionError((endpoint, op, v))


def states_matrix(frame, params):
    return np.column_stack([frame[e].map(lambda s: property_state(s, e, **params)).to_numpy(np.int8) for e in ENDPOINTS])


def joint(states, subset):
    sub = states[:, list(subset)]
    return (sub == 1).all(axis=1).astype(np.int8), (sub != 0).all(axis=1).astype(np.int8)


def bounds(a, b, lo, hi, k):
    ex_a, ex_b = a & ~b, b & ~a
    return ((lo[ex_a].sum() - hi[ex_b].sum()) / k,
            (hi[ex_a].sum() - lo[ex_b].sum()) / k,
            int(ex_a.sum()), int(ex_b.sum()))


def select_from_scores(values, names, fraction):
    order = np.lexsort((names, -values))
    k = int(math.ceil(fraction * len(names)))
    a = np.zeros(len(names), dtype=bool)
    a[order[:k]] = True
    return a, k


summary = {'generated_utc': datetime.now(timezone.utc).isoformat(),
           'inputs_sha256': {'later_raw': SHA(DATA / 'expansion_data_raw.csv'),
                             'later_scores': SHA(SCORES / 'policy_scores.csv'),
                             'dev_membership': SHA(DEV / 'selection_membership.csv'),
                             'train_adjudicated': SHA(TRAIN),
                             'structure_membership': SHA(STRUCT),
                             'script': SHA(Path(__file__))}}

# ------------------------------------------------------------------ later cohort
scores = pd.read_csv(SCORES / 'policy_scores.csv')
rules = [c for c in scores.columns if c not in ('Molecule Name', 'SMILES')]
assert len(rules) == 10
names = scores['Molecule Name'].to_numpy()
raw = pd.read_csv(DATA / 'expansion_data_raw.csv', dtype=str, keep_default_na=False)
raw = raw.set_index('Molecule Name').loc[names].reset_index()
assert raw['Molecule Name'].tolist() == list(names)
states = states_matrix(raw, PRIMARY)
lo3, hi3 = joint(states, (0, 1, 2))
assert int(lo3.sum()) == 283 and int((hi3 == 0).sum()) == 1374 and int((lo3 != hi3).sum()) == 625
selections = {f: {r: select_from_scores(scores[r].to_numpy(), names, f)[0] for r in rules} for f in BUDGETS}
K_LATER = {0.05: 115, 0.10: 229, 0.20: 457}

PROFILES = {'LogD': (0,), 'KSOL': (1,), 'HLM CLint': (2,), 'LogD + KSOL': (0, 1),
            'LogD + HLM CLint': (0, 2), 'KSOL + HLM CLint': (1, 2), 'All three': (0, 1, 2)}
PROFILE_ORDER = list(PROFILES)


def profile_map(states_matrix_, selections_, budgets, K):
    rows = []
    for pname, subset in PROFILES.items():
        lo, hi = joint(states_matrix_, subset)
        rows.append(dict(cohort_K=K, profile=pname, decidable=int((lo == hi).sum()),
                         known_successes=int(lo.sum()), unknown=int((lo != hi).sum())))
        for f in budgets:
            nid = 0
            for ra, rb in itertools.combinations(rules, 2):
                L, U, _, _ = bounds(selections_[f][ra], selections_[f][rb], lo, hi, K[f])
                nid += int(L > 1e-12 or U < -1e-12)
            rows[-1][f'identified_{int(round(f * 100))}pct'] = nid
    return pd.DataFrame(rows)


later_map = profile_map(states, selections, BUDGETS, K_LATER)
later_map.to_csv(OUT / 'profile_map_later.csv', index=False)

# ------------------------------------------------------------------ development cohort
adj = pd.read_csv(TRAIN)
dev_names = adj.loc[adj.subset == 'development', 'Molecule Name'].to_numpy()
train_raw = pd.read_csv(ROOT / 'artifacts/expansionrx_train_feasibility/raw_training_subset.csv',
                        dtype=str, keep_default_na=False)
train_raw = train_raw.set_index('Molecule Name').loc[dev_names].reset_index()
dev_states = states_matrix(train_raw, PRIMARY)
adj_dev = adj.loc[adj.subset == 'development'].set_index('Molecule Name').loc[dev_names]
assert (dev_states == adj_dev[['state_0', 'state_1', 'state_2']].to_numpy(np.int8)).all(), 'development state mismatch'
dev_lo, dev_hi = joint(dev_states, (0, 1, 2))
dev_members = pd.read_csv(DEV / 'selection_membership.csv')
dev_pos = {n: i for i, n in enumerate(dev_names)}
dev_sel = {}
for f in BUDGETS:
    dev_sel[f] = {}
    for r in rules:
        z = dev_members[(np.isclose(dev_members.fraction, f)) & (dev_members.rule == r) & (dev_members.full_selected)]
        a = np.zeros(len(dev_names), dtype=bool)
        a[[dev_pos[m] for m in z.molecule]] = True
        dev_sel[f][r] = a
    k = int(math.ceil(f * len(dev_names)))
    assert all(int(a.sum()) == k for a in dev_sel[f].values()), (f, [int(a.sum()) for a in dev_sel[f].values()])
K_DEV = {f: int(math.ceil(f * len(dev_names))) for f in BUDGETS}
dev_map = profile_map(dev_states, dev_sel, BUDGETS, K_DEV)
dev_map.to_csv(OUT / 'profile_map_development.csv', index=False)
summary['development'] = {'n': int(len(dev_names)), 'budgets': K_DEV,
                          'decidable': int((dev_lo == dev_hi).sum()), 'known_successes': int(dev_lo.sum()),
                          'unknown': int((dev_lo != dev_hi).sum()),
                          'hlm_missing': int((train_raw['HLM CLint'] == '').sum()),
                          'later_hlm_missing': int((raw['HLM CLint'] == '').sum())}

later_ident = later_map.set_index('profile')[['identified_5pct', 'identified_10pct', 'identified_20pct']].to_dict('index')
dev_ident = dev_map.set_index('profile')[['identified_5pct', 'identified_10pct', 'identified_20pct']].to_dict('index')
summary['profile_map'] = {'later': later_ident, 'development': dev_ident}
print('LATER MAP', json.dumps(later_ident))
print('DEV MAP', json.dumps(dev_ident))
print('DEV unknown', summary['development']['unknown'], 'HLM missing', summary['development']['hlm_missing'],
      'later HLM missing', summary['development']['later_hlm_missing'])

# ------------------------------------------------------------------ B coverage requirement + D witnesses
cov_rows, witness_rows = [], []
for f in BUDGETS:
    k = K_LATER[f]
    for ra, rb in itertools.combinations(rules, 2):
        a, b = selections[f][ra], selections[f][rb]
        L, U, n_ex_a, n_ex_b = bounds(a, b, lo3, hi3, k)
        ex_a, ex_b = a & ~b, b & ~a
        d = int(lo3[ex_a].sum() - lo3[ex_b].sum())
        u_a = int((lo3[ex_a] != hi3[ex_a]).sum())
        u_b = int((lo3[ex_b] != hi3[ex_b]).sum())
        # sharpness witnesses: completion attaining each endpoint. Base = known successes (1) and
        # known failures (0); only unresolved outcomes are assigned, and they are assigned to the
        # extremal values on the policy-exclusive sets.
        unknown = lo3 != hi3
        yL = lo3.astype(np.int8).copy()
        yL[ex_b & unknown] = 1                     # B-only unknowns: worst case for A - B
        yU = lo3.astype(np.int8).copy()
        yU[ex_a & unknown] = 1                     # A-only unknowns: best case for A - B
        achieved_L = (yL[a].sum() - yL[b].sum()) / k
        achieved_U = (yU[a].sum() - yU[b].sum()) / k
        assert abs(achieved_L - L) < 1e-12 and abs(achieved_U - U) < 1e-12, (ra, rb, f, achieved_L, L, achieved_U, U)
        for label, y, achieved, bound in [('lower', yL, achieved_L, L), ('upper', yU, achieved_U, U)]:
            witness_rows.append(dict(budget=f, policy_a=ra, policy_b=rb, extremum=label,
                                     bound=bound, achieved=achieved,
                                     assignment_sha256=hashlib.sha256(bytes(y)).hexdigest()))
        if abs(f - 0.10) > 1e-12:
            continue
        if d > 0:
            req, side, n_ex, u_opp, cur_known = max(0, u_b - d + 1), 'B', n_ex_b, u_b, (n_ex_b - u_b) / n_ex_b
        elif d < 0:
            req, side, n_ex, u_opp, cur_known = max(0, u_a + d + 1), 'A', n_ex_a, u_a, (n_ex_a - u_a) / n_ex_a
        else:
            req, side, n_ex, u_opp, cur_known = None, 'none', None, None, None
        cov_rows.append(dict(policy_a=ra, policy_b=rb, d=d, side=side,
                             opposite_exclusive=n_ex, opposite_unknown=u_opp,
                             required_resolved=req,
                             max_unknown_allowed=(abs(d) - 1) if d != 0 else None,
                             required_known_fraction=(1 - (abs(d) - 1) / n_ex) if d != 0 else None,
                             current_known_fraction=cur_known,
                             attainable=req is not None))
cov = pd.DataFrame(cov_rows)
cov.to_csv(OUT / 'coverage_requirement_primary45.csv', index=False)
wit = pd.DataFrame(witness_rows)
wit.to_csv(OUT / 'sharpness_witnesses.csv', index=False)
attain = cov[cov.attainable == True]
summary['coverage_requirement'] = {
    'pairs': int(len(cov)), 'attainable': int(len(attain)),
    'required_known_fraction_median': float(attain.required_known_fraction.median()),
    'required_known_fraction_min': float(attain.required_known_fraction.min()),
    'required_known_fraction_max': float(attain.required_known_fraction.max()),
    'current_known_fraction_median': float(attain.current_known_fraction.median()),
    'pairs_current_below_required': int((attain.current_known_fraction < attain.required_known_fraction).sum()),
    'required_resolved_median': float(attain.required_resolved.median()),
}
print('COVERAGE', json.dumps(summary['coverage_requirement']))
summary['sharpness'] = {'witness_rows': int(len(wit)), 'max_abs_deviation': float((wit.achieved - wit.bound).abs().max())}
print('WITNESSES', summary['sharpness'])

# ------------------------------------------------------------------ C measurement-value frontier
# universe: jointly unknown later compounds; each pair needs r_p of its opposite-side unknown records
unknown_idx = np.flatnonzero(lo3 != hi3)
universe = [names[i] for i in unknown_idx]
pos = {n: i for i, n in enumerate(universe)}
targets = {}
for row in attain.itertuples():
    a, b = selections[0.10][row.policy_a], selections[0.10][row.policy_b]
    ex = (b & ~a) if row.side == 'B' else (a & ~b)   # opposite-side exclusive set of the observed-favouring direction
    members = [pos[names[i]] for i in np.flatnonzero(ex & (lo3 != hi3))]
    targets[(row.policy_a, row.policy_b)] = (set(members), int(row.required_resolved))
pairs = list(targets)
chosen, covered, satisfied, curve = [], {p: 0 for p in pairs}, set(), []
remaining = {p: targets[p][1] for p in pairs}
while len(satisfied) < len(pairs):
    best, best_gain = None, -1.0
    for j in range(len(universe)):
        if j in chosen:
            continue
        gain = sum(1.0 / remaining[p] for p in pairs if p not in satisfied and j in targets[p][0])
        if gain > best_gain:
            best, best_gain = j, gain
    assert best is not None and best_gain > 0, (len(chosen), len(satisfied))
    chosen.append(best)
    for p in pairs:
        if p in satisfied:
            continue
        if best in targets[p][0]:
            covered[p] += 1
            remaining[p] = max(1, targets[p][1] - covered[p])
        if covered[p] >= targets[p][1]:
            satisfied.add(p)
    curve.append(dict(measurements=len(chosen), satisfied=len(satisfied)))
frontier = pd.DataFrame(curve)
frontier.to_csv(OUT / 'measurement_frontier.csv', index=False)
first = {}
for m in [1, 2, 5, 10, 20, 30, 40, 44]:
    z = frontier[frontier.satisfied >= m]
    first[m] = int(z.measurements.min()) if len(z) else None
summary['measurement_frontier'] = {'pairs_with_requirement': len(pairs), 'universe_unknown': len(universe),
                                   'measurements_for_m_pairs': first,
                                   'measurements_for_all': int(frontier.measurements.max()),
                                   'greedy_note': 'greedy plan; a plan of this size is achievable in the best case, so the true minimum is not larger.'}
print('FRONTIER', json.dumps(summary['measurement_frontier']))

# ------------------------------------------------------------------ E structural clustering of missingness
struct = pd.read_csv(STRUCT)
struct = struct[struct.official_split == 'test'].set_index('Molecule Name').loc[names].reset_index()
assert struct['Molecule Name'].tolist() == list(names)
hlm_missing = (raw['HLM CLint'] == '').to_numpy()
group = struct['murcko_nonstereo'].fillna('').to_numpy()
df = pd.DataFrame({'g': group, 'missing': hlm_missing, 'name': names})
df = df[df.g != '']
agg = df.groupby('g').agg(n=('missing', 'size'), miss=('missing', 'sum')).query('n >= 5')
agg['rate'] = agg.miss / agg.n
observed_share = float(df.merge(agg.reset_index(), on='g', suffixes=('', '_agg')).query('rate >= 0.8').missing.sum() / df.missing.sum()) if len(agg) else float('nan')
rng = np.random.RandomState(42)
perm = []
for _ in range(500):
    shuffled = df.copy()
    shuffled['missing'] = rng.permutation(shuffled['missing'].to_numpy())
    a2 = shuffled.groupby('g').agg(n=('missing', 'size'), miss=('missing', 'sum')).query('n >= 5')
    a2['rate'] = a2.miss / a2.n
    m2 = shuffled.merge(a2.reset_index(), on='g', suffixes=('', '_agg'))
    perm.append(float(m2.query('rate >= 0.8').missing.sum() / max(1, shuffled.missing.sum())))
perm = np.array(perm)
# leave-one-out group missing rate for selected compounds (primary budget, pooled over policies)
sel_any = np.zeros(len(names), dtype=bool)
for r in rules:
    sel_any |= selections[0.10][r]
gmap = df.set_index('name')
def loo_rate(n):
    if n not in gmap.index:
        return np.nan
    row = gmap.loc[n]
    if row.g not in agg.index:
        return np.nan
    grp = agg.loc[row.g]
    if grp.n <= 1:
        return np.nan
    return (grp.miss - int(row.missing)) / (grp.n - 1)
rates = np.array([loo_rate(n) for n in names])
pool_rate = float(np.nanmean(rates))
sel_rate = float(np.nanmean(rates[sel_any]))
struct_summary = {'scaffold_groups_ge5': int(len(agg)),
                  'groups_all_missing': int((agg.rate == 1).sum()),
                  'groups_none_missing': int((agg.rate == 0).sum()),
                  'share_missing_in_groups_rate_ge_0.8': observed_share,
                  'permutation_mean': float(perm.mean()), 'permutation_p95': float(np.percentile(perm, 95)),
                  'permutation_z': float((observed_share - perm.mean()) / (perm.std(ddof=1) + 1e-12)),
                  'pool_mean_loo_scaffold_missing_rate': pool_rate,
                  'selected_mean_loo_scaffold_missing_rate': sel_rate}
summary['missingness_structure'] = struct_summary
pd.DataFrame([struct_summary]).to_csv(OUT / 'missingness_structure_summary.csv', index=False)
print('STRUCTURE', json.dumps(struct_summary))

# ------------------------------------------------------------------ F naive complete-case magnitudes
pw = pd.read_csv(V4OUT / 'pairwise_identification.csv')
naive = []
for f in BUDGETS:
    z = pw[np.isclose(pw.fraction, f)]
    for col, name in [('restricted_delta', 'restricted_pool'), ('post_exact_delta', 'within_batch_exact')]:
        delta = z[col].dropna()
        pp = 100 * delta.abs()
        naive.append(dict(budget=f, report=name, non_tied=int((delta.abs() > 1e-12).sum()),
                          abs_delta_median_pp=float(pp[delta.abs() > 1e-12].median()),
                          abs_delta_p90_pp=float(pp[delta.abs() > 1e-12].quantile(0.9)),
                          abs_delta_max_pp=float(pp[delta.abs() > 1e-12].max())))
naive = pd.DataFrame(naive)
naive.to_csv(OUT / 'naive_report_magnitudes.csv', index=False)
print('NAIVE', naive.to_dict('records'))

summary['outputs'] = sorted(p.name for p in OUT.iterdir() if p.is_file())
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf-8')
print('DONE', OUT)
