import csv
import hashlib
import itertools
import json
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from statistics import median

out = Path(__file__).resolve().parent
revision = out.parent
root = revision.parent
source = root / 'artifacts/expansionrx_later_identification_v1'
read = lambda p: list(csv.DictReader(p.open(encoding='utf-8-sig', newline='')))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def write_csv(name, rows):
    with (revision / 'tables' / name).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

records = read(source / 'scored_outcomes.csv')
raw = {r['Molecule Name']: r for r in read(root / 'data/expansionrx_openadmet/expansion_data_raw.csv')}
later = read(root / 'data/expansionrx_openadmet/expansion_data_test.csv')
original = read(source / 'pairwise_identification.csv')
membership = read(source / 'selection_membership.csv')
endpoints = ['LogD', 'KSOL', 'HLM CLint']
policies = list(records[0])[2:12]
assert len(policies) == 10 and len(records) == 2282
assert [r['Molecule Name'] for r in records] == [r['Molecule Name'] for r in later]

def classify(text, j):
    text = text.strip()
    if not text:
        return -1
    op = text[0] if text[0] in '<>' else '='
    value = Decimal(text[1:].strip() if op != '=' else text)
    assert j == 0 or value >= 0
    if op == '=':
        return int(1 <= value <= 3) if j == 0 else int(value > 50) if j == 1 else int(value < 10)
    if j == 0:
        return 0 if (op == '>' and value >= 3) or (op == '<' and value <= 1) else -1
    if j == 1:
        return 1 if op == '>' and value >= 50 else 0 if op == '<' and value <= 50 else -1
    return 1 if op == '<' and 0 < value <= 10 else 0 if op == '>' and value >= 10 else -1

labels, states, exact = {}, {}, {}
for row, test in zip(records, later):
    identifier = row['Molecule Name']
    assert row['SMILES'] == raw[identifier]['SMILES'] == test['SMILES']
    for endpoint in endpoints:
        assert row['raw_' + endpoint] == raw[identifier][endpoint]
    ss = tuple(classify(row['raw_' + e], j) for j, e in enumerate(endpoints))
    assert ss == tuple(int(row[f'state_{j}']) for j in range(3))
    lo, hi = int(all(s == 1 for s in ss)), int(all(s != 0 for s in ss))
    assert (lo, hi) == (int(float(row['success_lower'])), int(float(row['success_upper'])))
    ee = all(row['raw_' + e].strip() and row['raw_' + e].strip()[0] not in '<>' for e in endpoints)
    assert ee == (row['all_exact'] == 'True')
    labels[identifier], states[identifier], exact[identifier] = (lo, hi), ss, ee

unknown = {i for i, v in labels.items() if v == (0, 1)}
patterns = Counter(tuple(s == -1 for s in states[i]) for i in unknown)
pattern_rows = []
for pattern, n in sorted(patterns.items(), key=lambda item: (-item[1], item[0])):
    pattern_rows.append(dict(unknown_endpoints=' + '.join(e for e, u in zip(endpoints, pattern) if u),
                             LogD_unknown=pattern[0], KSOL_unknown=pattern[1], HLM_unknown=pattern[2], count=n,
                             fraction_of_joint_unknown=n / len(unknown)))
write_csv('joint_unknown_endpoint_patterns.csv', pattern_rows)

pairs = []
selection_checks = 0
for q in [Fraction(1, 20), Fraction(1, 10), Fraction(1, 5)]:
    k = (q.numerator * len(records) + q.denominator - 1) // q.denominator
    selected = {p: {r['Molecule Name'] for r in sorted(records, key=lambda r: (-Decimal(r[p]), r['Molecule Name']))[:k]} for p in policies}
    for policy in policies:
        saved = {r['molecule'] for r in membership if Fraction(r['fraction']) == q and r['rule'] == policy and r['full_selected'] == 'True'}
        assert selected[policy] == saved
        selection_checks += 1
    for a, b in itertools.combinations(policies, 2):
        shared, only_a, only_b = selected[a] & selected[b], selected[a] - selected[b], selected[b] - selected[a]
        sa, sb = sum(labels[i][0] for i in only_a), sum(labels[i][0] for i in only_b)
        ua, ub = len(only_a & unknown), len(only_b & unknown)
        d = sa - sb
        lo, hi = Fraction(d - ub, k), Fraction(d + ua, k)
        saved = next(r for r in original if Fraction(r['fraction']) == q and r['policy_a'] == a and r['policy_b'] == b)
        assert abs(float(lo) - float(saved['lower'])) < 1e-12
        assert abs(float(hi) - float(saved['upper'])) < 1e-12
        assert len(shared) == int(saved['overlap']) and ua + ub == int(saved['unique_unknown'])
        assert hi - lo == Fraction(ua + ub, k) and lo < 0 < hi
        pairs.append(dict(fraction=float(q), k=k, policy_a=a, policy_b=b, overlap=len(shared),
                          shared_unknown=len(shared & unknown), exclusive_per_side=len(only_a),
                          known_success_a_only=sa, known_success_b_only=sb, known_success_difference=d,
                          unknown_a_only=ua, unknown_b_only=ub, exclusive_unknown_total=ua + ub,
                          lower_numerator=d - ub, upper_numerator=d + ua, lower_fraction=str(lo), upper_fraction=str(hi),
                          lower_pp=100 * float(lo), upper_pp=100 * float(hi), width_pp=100 * float(hi - lo)))
write_csv('paired_decomposition_all_135.csv', pairs)
primary = [r for r in pairs if r['fraction'] == .1]
write_csv('paired_decomposition_primary_45.csv', primary)
ordered = sorted(primary, key=lambda r: (r['exclusive_unknown_total'], r['policy_a'], r['policy_b']))
representative = [dict(selection=label, width_rank=rank + 1, **ordered[rank]) for label, rank in [('Minimum width', 0), ('Median-ranked width', 22), ('Maximum width', 44)]]
write_csv('representative_pairs.csv', representative)
summary = dict(scope='Post-review exploratory decomposition and independent saved-result recomputation; no retraining or changed outcomes.',
               n=len(records), complete_exact=sum(exact.values()), joint_decidable=len(records) - len(unknown),
               known_success=sum(x[0] for x in labels.values()), known_failure=sum(x == (0, 0) for x in labels.values()),
               joint_unknown=len(unknown), successes_with_valid_inequalities=sum(labels[i][0] and not exact[i] for i in labels),
               joint_unknown_involving_HLM=sum(states[i][2] == -1 for i in unknown),
               joint_unknown_only_HLM=sum(states[i] == (1, 1, -1) for i in unknown),
               endpoint_states_matched=3 * len(records), saved_selection_sets_matched=selection_checks,
               pair_intervals_matched=len(pairs), all_135_strictly_cross_zero=True,
               primary=dict(k=229, count=len(primary),
                            width_pp_min_median_max=[min(r['width_pp'] for r in primary), median(r['width_pp'] for r in primary), max(r['width_pp'] for r in primary)],
                            overlap_min_median_max=[min(r['overlap'] for r in primary), median(r['overlap'] for r in primary), max(r['overlap'] for r in primary)],
                            exclusive_unknown_min_median_max=[min(r['exclusive_unknown_total'] for r in primary), median(r['exclusive_unknown_total'] for r in primary), max(r['exclusive_unknown_total'] for r in primary)]),
               representative_selection='Sort 45 primary pairs by exclusive unknown count, then policy_a and policy_b lexicographically; take ranks 1, 23, 45. Orientation remains saved policy_a minus policy_b.',
               representative_pairs=representative,
               source_sha256={str(p.relative_to(root)): sha(p) for p in [source / 'scored_outcomes.csv', source / 'pairwise_identification.csv', source / 'selection_membership.csv', root / 'data/expansionrx_openadmet/expansion_data_raw.csv', root / 'data/expansionrx_openadmet/expansion_data_test.csv']})
assert [summary[key] for key in ['complete_exact', 'joint_decidable', 'known_success', 'known_failure', 'joint_unknown', 'successes_with_valid_inequalities', 'joint_unknown_involving_HLM', 'joint_unknown_only_HLM']] == [776, 1657, 283, 1374, 625, 153, 623, 567]
(out / 'validation.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
