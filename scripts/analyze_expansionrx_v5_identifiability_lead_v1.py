"""Identifiability analysis using saved scores, selections and reported measurements.

Analyses:
  1. Check the summary results.
  2. Endpoint-profile identifiability map (7 non-empty subsets of LogD/KSOL/HLM x 3 budgets).
  3. Pre-specified threshold sensitivity grid (10 distinct scenarios).
  4. Identification distance for the 45 primary-budget pairs.
"""
import hashlib
import json
import math
import re
import itertools
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/expansionrx_openadmet'
SCORES = ROOT / 'artifacts/expansionrx_later_fixed_scores_v2'
V4OUT = ROOT / 'artifacts/expansionrx_later_identification_v1'
OUT = ROOT / 'artifacts/expansionrx_v5_identifiability_lead_v1'
OUT.mkdir(exist_ok=False)

ENDPOINTS = ['LogD', 'KSOL', 'HLM CLint']
BUDGETS = [0.05, 0.10, 0.20]
K = {0.05: 115, 0.10: 229, 0.20: 457}
PRIMARY = dict(logd_lo=1.0, logd_hi=3.0, ksol_min=50.0, hlm_max=10.0)
SHA = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
NUM_RE = re.compile(r'\s*(<=|>=|<|>|=)?\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*')


def property_state(value, endpoint, logd_lo, logd_hi, ksol_min, hlm_max):
    """Classify an endpoint measurement: -1 unknown, 0 fail, 1 pass."""
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


def states_matrix(raw, params):
    cols = [raw[e].map(lambda s: property_state(s, e, **params)).to_numpy(np.int8) for e in ENDPOINTS]
    return np.column_stack(cols)


def joint_bounds(states, subset):
    """Per-compound [l,u] for a profile = subset of endpoint indices."""
    sub = states[:, list(subset)]
    lo = (sub == 1).all(axis=1).astype(np.int8)
    hi = (sub != 0).all(axis=1).astype(np.int8)
    return lo, hi


def select(rule_scores, names, fraction):
    order = np.lexsort((names, -rule_scores))
    k = int(math.ceil(fraction * len(names)))
    a = np.zeros(len(names), dtype=bool)
    a[order[:k]] = True
    return a, k


def pair_bounds(a, b, lo, hi, k):
    ex_a = a & ~b
    ex_b = b & ~a
    lower = (lo[ex_a].sum() - hi[ex_b].sum()) / k
    upper = (hi[ex_a].sum() - lo[ex_b].sum()) / k
    return float(lower), float(upper), int(ex_a.sum()), int(ex_b.sum())


def identified(lower, upper, tol=1e-12):
    return bool(lower > tol or upper < -tol)


def touch_zero(lower, upper, tol=1e-12):
    return bool(abs(lower) < tol or abs(upper) < tol)


# ---------------------------------------------------------------- load inputs
scores = pd.read_csv(SCORES / 'policy_scores.csv')
rules = [c for c in scores.columns if c not in ('Molecule Name', 'SMILES')]
assert len(rules) == 10, rules
names = scores['Molecule Name'].to_numpy()
raw = pd.read_csv(DATA / 'expansion_data_raw.csv', dtype=str, keep_default_na=False)
raw = raw.set_index('Molecule Name').loc[names].reset_index()
assert raw['Molecule Name'].tolist() == list(names)
assert scores['SMILES'].tolist() == raw['SMILES'].tolist()

# saved selections must reproduce from scores
rank = pd.read_csv(SCORES / 'rankings_and_selection.csv')
rank10 = rank[rank.top10]
for rule in rules:
    a, _ = select(scores[rule].to_numpy(), names, 0.10)
    saved = set(rank10.loc[rank10.policy == rule, 'Molecule Name'])
    assert set(names[a]) == saved, rule
selections = {f: {r: select(scores[r].to_numpy(), names, f)[0] for r in rules} for f in BUDGETS}

summary = {'generated_utc': datetime.now(timezone.utc).isoformat(), 'n_later': int(len(scores)),
           'input_sha256': {'policy_scores': SHA(SCORES / 'policy_scores.csv'),
                            'rankings_and_selection': SHA(SCORES / 'rankings_and_selection.csv'),
                            'raw_strings': SHA(DATA / 'expansion_data_raw.csv'),
                            'script_sha256': SHA(Path(__file__))}}

# ---------------------------------------------------------------- 1 reproduction gate
states = states_matrix(raw, PRIMARY)
lo3, hi3 = joint_bounds(states, (0, 1, 2))
# v4 complete-exact: finite numeric for all three endpoints, nonnegative for KSOL/HLM, decidable states
numeric = np.column_stack([pd.to_numeric(raw[e], errors='coerce').to_numpy(float) for e in ENDPOINTS])
exact = np.isfinite(numeric).all(axis=1) & np.all(numeric[:, 1:] >= 0, axis=1) & (states >= 0).all(axis=1)

succ = int(lo3.sum())
fail = int((hi3 == 0).sum())
unk = int((lo3 != hi3).sum())
decidable = int((lo3 == hi3).sum())
succ_ineq = int(((lo3 == 1) & ~exact).sum())
raw_nonnull = raw[ENDPOINTS].ne('').to_numpy()
cov = {}
for j, e in enumerate(ENDPOINTS):
    col = raw[e]
    exact_n = int(pd.to_numeric(col.replace('', np.nan), errors='coerce').notna().sum())
    ineq_n = int(col.str.match(r'^\s*(<=|>=|<|>)\s*[0-9.]').sum())
    missing_n = int((col == '').sum())
    semantic_invalid = 0  # asserted separately in v4; observed zero in later cohort
    cov[e] = dict(exact=exact_n, inequality=ineq_n, missing=missing_n, semantic_invalid=semantic_invalid, total=int(len(col)))

gate = {
    'n': 2282, 'complete_exact': int(exact.sum()), 'decidable': decidable,
    'known_successes': succ, 'known_failures': fail, 'unknown': unk,
    'successes_with_inequality': succ_ineq, 'successes_complete_exact': int((lo3[exact] == 1).sum()),
    'endpoint_coverage': cov,
}
assert gate['n'] == 2282
assert gate['complete_exact'] == 776, gate['complete_exact']
assert gate['decidable'] == 1657, gate['decidable']
assert gate['known_successes'] == 283, gate['known_successes']
assert gate['known_failures'] == 1374, gate['known_failures']
assert gate['unknown'] == 625, gate['unknown']
assert gate['successes_with_inequality'] == 153, gate['successes_with_inequality']
assert gate['successes_complete_exact'] == 130, gate['successes_complete_exact']
assert cov['LogD']['exact'] == 2270 and cov['LogD']['inequality'] == 2 and cov['LogD']['missing'] == 10, cov['LogD']
assert cov['KSOL']['exact'] == 2170 and cov['KSOL']['inequality'] == 11 and cov['KSOL']['missing'] == 101, cov['KSOL']
assert cov['HLM CLint']['exact'] == 782 and cov['HLM CLint']['inequality'] == 228 and cov['HLM CLint']['missing'] == 1272, cov['HLM CLint']

# 1b. what actually creates decidability beyond complete-exact panels:
# counterfactual with every inequality-derived state set to unknown.
INQ = re.compile(r'^\s*(<=|>=|<|>)\s*')
raw_exact_only = raw.copy()
for e in ENDPOINTS:
    mask = raw_exact_only[e].str.match(INQ)
    raw_exact_only.loc[mask, e] = ''
states_no_ineq = states_matrix(raw_exact_only, PRIMARY)
lo_no, hi_no = joint_bounds(states_no_ineq, (0, 1, 2))
dec_full = lo3 == hi3
dec_no = lo_no == hi_no
ineq_dependent = dec_full & ~dec_no
gate['decidable_full'] = int(dec_full.sum())
gate['decidable_without_inequalities'] = int(dec_no.sum())
gate['decidable_dependent_on_inequalities'] = int(ineq_dependent.sum())
gate['ineq_dependent_successes'] = int((ineq_dependent & (lo3 == 1)).sum())
gate['ineq_dependent_failures'] = int((ineq_dependent & (hi3 == 0)).sum())
# records decidable but not complete-exact, and not inequality dependent
gate['nonexact_decidable_not_ineq'] = int((dec_full & ~exact & ~ineq_dependent).sum())
assert gate['ineq_dependent_successes'] == 153, gate['ineq_dependent_successes']
assert gate['decidable_full'] == 1657
assert gate['decidable_full'] == gate['decidable_without_inequalities'] + gate['decidable_dependent_on_inequalities']

# 1c. selection enriches unknown outcomes: unknown rate inside each selected batch vs pool.
pool_unknown_rate = unk / len(scores)
sel_rows = []
for f in BUDGETS:
    for rule in rules:
        a = selections[f][rule]
        u = int((lo3[a] != hi3[a]).sum())
        sel_rows.append(dict(budget=f, k=K[f], policy=rule, selected=int(a.sum()), unknown=u,
                             unknown_rate=u / K[f], pool_unknown_rate=pool_unknown_rate,
                             enrichment=(u / K[f]) / pool_unknown_rate))
gate['selected_unknown_rate_min'] = float(min(r['unknown_rate'] for r in sel_rows))
gate['selected_unknown_rate_max'] = float(max(r['unknown_rate'] for r in sel_rows))
gate['selected_unknown_rate_mean'] = float(np.mean([r['unknown_rate'] for r in sel_rows]))
gate['pool_unknown_rate'] = float(pool_unknown_rate)
pd.DataFrame(sel_rows).to_csv(OUT / 'selected_unknown_rates.csv', index=False)

# 1d. Table 1 source: mutually exclusive decidability classes and endpoint coverage.
cls = pd.DataFrame([
    dict(cls='Complete-exact panel', records=int(exact.sum()),
         known_successes=int((lo3[exact] == 1).sum()), known_failures=int((hi3[exact] == 0).sum())),
    dict(cls='Decidable through an exact failure, other endpoints missing',
         records=int(gate['decidable_full'] - exact.sum() - gate['decidable_dependent_on_inequalities']),
         known_successes=0,
         known_failures=int(gate['decidable_full'] - exact.sum() - gate['decidable_dependent_on_inequalities'])),
    dict(cls='Decidable only through a valid inequality', records=int(gate['decidable_dependent_on_inequalities']),
         known_successes=int(gate['ineq_dependent_successes']), known_failures=int(gate['ineq_dependent_failures'])),
    dict(cls='Joint outcome unknown', records=int(unk), known_successes=0, known_failures=0),
])
assert cls.records.sum() == 2282 and cls.known_successes.sum() == 283 and cls.known_failures.sum() == 1374
assert cls.loc[1, 'records'] == 720 and cls.loc[2, 'records'] == 161
cls.to_csv(OUT / 'table_decidability_decomposition.csv', index=False)
covdf = pd.DataFrame([dict(endpoint=e, **v) for e, v in cov.items()])
covdf.to_csv(OUT / 'table_endpoint_coverage.csv', index=False)

# primary pairwise results (full 3-endpoint profile)
prim = []
for f in BUDGETS:
    for ra, rb in itertools.combinations(rules, 2):
        a, b = selections[f][ra], selections[f][rb]
        low, up, _, _ = pair_bounds(a, b, lo3, hi3, K[f])
        prim.append(dict(budget=f, k=K[f], policy_a=ra, policy_b=rb, lower=low, upper=up,
                         identified=identified(low, up), touches_zero=touch_zero(low, up)))
prim = pd.DataFrame(prim)
gate['primary_pairs_identified'] = int(prim.identified.sum())
gate['primary_pairs_touching_zero'] = int(prim.touches_zero.sum())
assert int(prim.identified.sum()) == 0, prim[prim.identified]
w10 = (prim.budget == 0.10)
widths_pp = 100 * (prim.loc[w10, 'upper'] - prim.loc[w10, 'lower'])
gate['primary_width_min_pp'] = round(float(widths_pp.min()), 6)
gate['primary_width_median_pp'] = round(float(widths_pp.median()), 6)
gate['primary_width_max_pp'] = round(float(widths_pp.max()), 6)
assert abs(gate['primary_width_min_pp'] - 5.240174672489083) < 1e-9 or abs(gate['primary_width_min_pp'] - 5.24) < 0.01
assert abs(gate['primary_width_max_pp'] - 62.00873362445415) < 1e-9 or abs(gate['primary_width_max_pp'] - 62.01) < 0.01
# width == exclusive unknown / k identity
for f in BUDGETS:
    for ra, rb in itertools.combinations(rules, 2):
        a, b = selections[f][ra], selections[f][rb]
        ex = a ^ b
        u = int(((lo3 != hi3) & ex).sum())
        row = prim[(prim.budget == f) & (prim.policy_a == ra) & (prim.policy_b == rb)].iloc[0]
        assert abs((row.upper - row.lower) - u / K[f]) < 1e-12
prim.to_csv(OUT / 'reproduction_primary_pairs.csv', index=False)
(OUT / 'reproduction_gate.json').write_text(json.dumps(gate, indent=2), encoding='utf-8')
print('REPRODUCTION GATE PASS', json.dumps({k: gate[k] for k in ('complete_exact', 'decidable', 'known_successes', 'unknown', 'successes_with_inequality', 'primary_pairs_identified', 'primary_width_min_pp', 'primary_width_median_pp', 'primary_width_max_pp')}))

# ---------------------------------------------------------------- 2 endpoint-profile map
profiles = {
    'LogD': (0,), 'KSOL': (1,), 'HLM': (2,),
    'LogD+KSOL': (0, 1), 'LogD+HLM': (0, 2), 'KSOL+HLM': (1, 2),
    'LogD+KSOL+HLM': (0, 1, 2),
}
profile_order = ['LogD', 'KSOL', 'HLM', 'LogD+KSOL', 'LogD+HLM', 'KSOL+HLM', 'LogD+KSOL+HLM']
rows_map, rows_pairs = [], []
for pname in profile_order:
    sub = profiles[pname]
    lo, hi = joint_bounds(states, sub)
    rows_map.append(dict(profile=pname, n_endpoints=len(sub), decidable=int((lo == hi).sum()),
                         known_successes=int(lo.sum()), known_failures=int((hi == 0).sum()),
                         unknown=int((lo != hi).sum())))
    for f in BUDGETS:
        for ra, rb in itertools.combinations(rules, 2):
            a, b = selections[f][ra], selections[f][rb]
            low, up, na, nb = pair_bounds(a, b, lo, hi, K[f])
            rows_pairs.append(dict(profile=pname, budget=f, k=K[f], policy_a=ra, policy_b=rb,
                                   lower=low, upper=up, identified=identified(low, up),
                                   width_pp=100 * (up - low), exclusive_a=na, exclusive_b=nb))
prof_pairs = pd.DataFrame(rows_pairs)
prof_map = pd.DataFrame(rows_map)
agg = []
for pname in profile_order:
    for f in BUDGETS:
        z = prof_pairs[(prof_pairs.profile == pname) & (prof_pairs.budget == f)]
        agg.append(dict(profile=pname, budget=f, k=K[f], identified=int(z.identified.sum()),
                        tied=int((~z.identified & ((z.lower == 0) & (z.upper == 0))).sum()),
                        unresolved=int((~z.identified).sum()),
                        width_min_pp=float(z.width_pp.min()), width_median_pp=float(z.width_pp.median()),
                        width_max_pp=float(z.width_pp.max())))
prof_agg = pd.DataFrame(agg)
prof_agg.to_csv(OUT / 'profile_identifiability.csv', index=False)
prof_pairs.to_csv(OUT / 'profile_pair_intervals.csv', index=False)
prof_map.to_csv(OUT / 'profile_outcome_counts.csv', index=False)
matrix = {p: {str(f): int(prof_agg[(prof_agg.profile == p) & (prof_agg.budget == f)].identified.iloc[0]) for f in BUDGETS} for p in profile_order}
summary['profile_map'] = {'identified_matrix': matrix, 'outcome_counts': prof_map.to_dict('records')}
print('PROFILE MAP', json.dumps(matrix))

# ---------------------------------------------------------------- 3 threshold sensitivity
GRID = [
    (1, 'primary', [1, 3], 50, 10),
    (2, 'logd', [0, 2], 50, 10),
    (3, 'logd', [2, 4], 50, 10),
    (4, 'logd', [1, 4], 50, 10),
    (5, 'ksol', [1, 3], 30, 10),
    (6, 'ksol', [1, 3], 100, 10),
    (7, 'hlm', [1, 3], 50, 5),
    (8, 'hlm', [1, 3], 50, 20),
    (9, 'loose_all', [0, 4], 30, 20),
    (10, 'strict_all', [1, 2], 100, 5),
]
grid_rows, grid_pairs = [], []
for sid, axis, window, kmin_, hmax in GRID:
    params = dict(logd_lo=float(window[0]), logd_hi=float(window[1]), ksol_min=float(kmin_), hlm_max=float(hmax))
    st = states_matrix(raw, params)
    lo, hi = joint_bounds(st, (0, 1, 2))
    grid_rows.append(dict(scenario=sid, axis=axis, logd_lo=window[0], logd_hi=window[1], ksol_min=kmin_, hlm_max=hmax,
                          decidable=int((lo == hi).sum()), known_successes=int(lo.sum()),
                          known_failures=int((hi == 0).sum()), unknown=int((lo != hi).sum())))
    for f in BUDGETS:
        nid = 0
        for ra, rb in itertools.combinations(rules, 2):
            a, b = selections[f][ra], selections[f][rb]
            low, up, _, _ = pair_bounds(a, b, lo, hi, K[f])
            ident = identified(low, up)
            nid += int(ident)
            grid_pairs.append(dict(scenario=sid, axis=axis, budget=f, k=K[f], policy_a=ra, policy_b=rb,
                                   lower=low, upper=up, identified=ident, width_pp=100 * (up - low)))
        grid_rows[-1][f'identified_{int(round(f * 100))}pct'] = nid
grid_df = pd.DataFrame(grid_rows)
grid_pairs = pd.DataFrame(grid_pairs)
grid_df.to_csv(OUT / 'threshold_sensitivity.csv', index=False)
grid_pairs.to_csv(OUT / 'threshold_pair_intervals.csv', index=False)
summary['threshold_sensitivity'] = {'grid': GRID, 'rows': grid_df.to_dict('records')}
print('THRESHOLD GRID identified@10%:', grid_df.set_index('scenario')['identified_10pct'].to_dict())

# ---------------------------------------------------------------- 4 identification distance
dist_rows = []
f = 0.10
for ra, rb in itertools.combinations(rules, 2):
    a, b = selections[f][ra], selections[f][rb]
    ex_a, ex_b = a & ~b, b & ~a
    d = int(lo3[ex_a].sum() - lo3[ex_b].sum())
    u_a = int((lo3[ex_a] != hi3[ex_a]).sum())
    u_b = int((lo3[ex_b] != hi3[ex_b]).sum())
    low, up, _, _ = pair_bounds(a, b, lo3, hi3, K[f])
    # resolving m of the B-only unknowns to failure raises lower by m/k (d fixed);
    # resolving m of the A-only unknowns to failure lowers upper by m/k.
    m_for_positive = max(0, u_b - d + 1) if d >= 1 else None
    m_for_negative = max(0, u_a + d + 1) if d <= -1 else None
    if m_for_positive is not None and m_for_positive > u_b:
        m_for_positive = None
    if m_for_negative is not None and m_for_negative > u_a:
        m_for_negative = None
    fav = 'A>B' if d > 0 else ('B>A' if d < 0 else 'tie')
    m_fav = m_for_positive if d > 0 else (m_for_negative if d < 0 else None)
    # how many of the exclusive unknowns are HLM-only unknown (LogD+KSOL pass, HLM unknown)
    hlm_only = (states[:, 0] == 1) & (states[:, 1] == 1) & (states[:, 2] == -1)
    dist_rows.append(dict(policy_a=ra, policy_b=rb, overlap=int((a & b).sum()), d=d, u_a=u_a, u_b=u_b,
                          lower=low, upper=up, observed_favouring=fav,
                          m_min_favouring=m_fav, attainable=bool(m_fav is not None),
                          exclusive_unknown=int(((lo3 != hi3) & (a ^ b)).sum()),
                          exclusive_unknown_HLM_only=int((hlm_only & (a ^ b)).sum())))
dist = pd.DataFrame(dist_rows)
dist.to_csv(OUT / 'identification_distance_primary45.csv', index=False)
summary['identification_distance'] = {
    'pairs': int(len(dist)), 'pairs_with_nonzero_d': int((dist.d != 0).sum()),
    'pairs_attainable_observed_favouring': int(dist.attainable.sum()),
    'm_min_median': float(dist.loc[dist.attainable, 'm_min_favouring'].median()),
    'm_min_min': int(dist.loc[dist.attainable, 'm_min_favouring'].min()),
    'm_min_max': int(dist.loc[dist.attainable, 'm_min_favouring'].max()),
    'exclusive_unknown_total': int(dist.exclusive_unknown.sum()),
    'exclusive_unknown_HLM_only_total': int(dist.exclusive_unknown_HLM_only.sum()),
}
print('IDENTIFICATION DISTANCE', json.dumps(summary['identification_distance']))

summary['outputs'] = sorted(p.name for p in OUT.iterdir() if p.is_file())
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False), encoding='utf-8')
print('DONE', OUT)
