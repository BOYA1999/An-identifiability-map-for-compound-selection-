"""v5 figures: identifiability map (Figure 2) and evidence gap (Figure 3).

Sources (all in tables/): profile_identifiability.csv, profile_outcome_counts.csv,
profile_pair_intervals.csv, joint_unknown_endpoint_patterns.csv,
identification_distance_primary45.csv, selected_unknown_rates.csv.
"""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

root = Path(__file__).resolve().parent.parent
figures = root / 'figures'
read = lambda name: list(csv.DictReader((root / 'tables' / name).open(encoding='utf-8')))

profile_agg = read('profile_identifiability.csv')
profile_outcomes = read('profile_outcome_counts.csv')
profile_pairs = read('profile_pair_intervals.csv')
patterns = read('joint_unknown_endpoint_patterns.csv')
distance = read('identification_distance_primary45.csv')
exposure = read('selected_unknown_rates.csv')

order = ['LogD', 'KSOL', 'HLM', 'LogD+KSOL', 'LogD+HLM', 'KSOL+HLM', 'LogD+KSOL+HLM']
short = {'LogD': 'LogD', 'KSOL': 'KSOL', 'HLM': 'HLM CLint', 'LogD+KSOL': 'LogD + KSOL',
         'LogD+HLM': 'LogD + HLM', 'KSOL+HLM': 'KSOL + HLM', 'LogD+KSOL+HLM': 'All three'}
budgets = [0.05, 0.10, 0.20]
colors = ['#3D6A80', '#A17453', '#697D68']
ink, muted = '#253744', '#66727A'
plt.style.use(Path(__file__).parent / 'academic.mplstyle')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8, 'axes.grid': False, 'axes.edgecolor': '#BBC3C8', 'text.color': ink,
                     'axes.labelcolor': ink, 'figure.dpi': 150, 'savefig.dpi': 300, 'svg.fonttype': 'none'})


def export(fig, stem, sources):
    for extension in ['png', 'pdf', 'svg']:
        options = {'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}
        fig.savefig(figures / f'{stem}.{extension}', dpi=300, bbox_inches='tight', pad_inches=.06, **options)
    plt.close(fig)
    return {'figure': stem, 'sources': sources,
            'exports': [f'figures/{stem}.{ext}' for ext in ['png', 'pdf', 'svg']],
            'source_sha256': {s: hashlib.sha256((root / 'tables' / s).read_bytes()).hexdigest() for s in sources}}


catalog_path = figures / 'figure_catalog.json'
catalog = json.loads(catalog_path.read_text()) if catalog_path.exists() else {}

# ---------------------------------------------------------------- Figure 2
matrix = np.zeros((len(order), len(budgets)))
for i, p in enumerate(order):
    for j, f in enumerate(budgets):
        row = [r for r in profile_agg if r['profile'] == p and abs(float(r['budget']) - f) < 1e-9][0]
        matrix[i, j] = int(row['identified'])

fig = plt.figure(figsize=(7.6, 4.9))
ax = fig.add_axes([.155, .40, .245, .48])
cmap = LinearSegmentedColormap.from_list('idmap', ['#F2F4F5', '#C6D3DA', '#7D9AAA', '#3D6A80'])
im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=45, aspect='auto')
for i in range(len(order)):
    for j in range(len(budgets)):
        value = int(matrix[i, j])
        ax.text(j, i, f'{value}', ha='center', va='center', fontsize=10,
                color='white' if value > 22 else ink, weight='bold')
ax.set_xticks(range(3), ['5%', '10%', '20%'])
ax.set_yticks(range(len(order)), [short[p] for p in order])
ax.set_xlabel('Selection budget')
ax.tick_params(axis='y', length=0, pad=4)
ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
ax.grid(which='minor', color='white', linewidth=1.6)
ax.tick_params(which='minor', length=0)
cbar = fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
cbar.set_ticks([0, 15, 30, 45])
cbar.set_label('Identified, of 45', fontsize=8)
cbar.ax.tick_params(labelsize=7)

axb = fig.add_axes([.50, .40, .21, .48])
counts = [int([r for r in profile_outcomes if r['profile'] == p][0]['unknown']) for p in order]
bars = axb.barh(range(len(order)), counts, color=['#A17453' if 'HLM' in p else '#879AA6' for p in order], height=.62)
for i, value in enumerate(counts):
    axb.text(value * 1.12, i, f'{value:,}', va='center', fontsize=7.5)
axb.set_xscale('log')
axb.set_xlim(7, 3000)
axb.set_xticks([10, 100, 1000], ['10', '100', '1,000'])
axb.set_yticks(range(len(order)), [])
axb.invert_yaxis()
axb.set_xlabel('Joint outcomes unknown (log scale)')
axb.spines['left'].set_visible(False)
axb.tick_params(axis='y', length=0)
axb.axhspan(-.5, 2.5, color='#3D6A80', alpha=.05, zorder=0)
axb.axhspan(3.5, 6.5, color='#A17453', alpha=.06, zorder=0)

axc = fig.add_axes([.80, .40, .185, .48])
data = []
for p in order:
    widths = [100 * (float(r['upper']) - float(r['lower'])) for r in profile_pairs
              if r['profile'] == p and abs(float(r['budget']) - 0.10) < 1e-9]
    data.append(widths)
box = axc.boxplot(data, vert=False, positions=range(len(order)), widths=.6, patch_artist=True, showfliers=False,
                  medianprops=dict(color='white', linewidth=1.1), whiskerprops=dict(color=muted, linewidth=.8),
                  capprops=dict(color=muted, linewidth=.8), boxprops=dict(linewidth=.6))
for patch, p in zip(box['boxes'], order):
    patch.set_facecolor('#A17453' if 'HLM' in p else '#879AA6')
    patch.set_alpha(.85)
axc.set_yticks(range(len(order)), [''] * len(order))
axc.invert_yaxis()
axc.set_ylim(len(order) - .5, -.5)
axc.set_xlim(0, 85)
axc.set_xticks([0, 40, 80])
axc.set_xlabel('Interval width, pp')
axc.grid(axis='x', color='#E4E8EB', linewidth=.5)

fig.text(.02, .965, 'a', fontsize=12, weight='bold')
fig.text(.055, .965, 'Which comparisons are identified?', fontsize=10, weight='bold')
fig.text(.50, .965, 'b', fontsize=12, weight='bold')
fig.text(.532, .965, 'Unknown outcomes', fontsize=10, weight='bold')
fig.text(.80, .965, 'c', fontsize=12, weight='bold')
fig.text(.832, .965, 'Width at 10%', fontsize=10, weight='bold')
fig.text(.5, .045, 'Every profile containing HLM CLint identifies none; profiles without it identify most.\n'
                   'Profiles are different scenarios, not alternative estimates of one rate.', ha='center', fontsize=8, color=muted, linespacing=1.3)
catalog['figure_2_identifiability_map'] = export(
    fig, 'figure_2_identifiability_map',
    ['profile_identifiability.csv', 'profile_outcome_counts.csv', 'profile_pair_intervals.csv'])
catalog['figure_2_identifiability_map']['claim'] = (
    'Within this fixed pool, every HLM-inclusive profile identifies 0 of 45 comparisons, '
    'every HLM-excluding profile 38-45.')

# ---------------------------------------------------------------- Figure 3
fig = plt.figure(figsize=(7.6, 3.8))
left = fig.add_axes([.17, .23, .25, .58])
pattern_short = {'HLM CLint': 'HLM only', 'KSOL + HLM CLint': 'KSOL + HLM', 'LogD + HLM CLint': 'LogD + HLM',
                 'LogD + KSOL + HLM CLint': 'All three', 'LogD + KSOL': 'LogD + KSOL', 'KSOL': 'KSOL only'}
left.barh(range(6), [int(r['count']) for r in patterns],
          color=['#A17453' if r['HLM_unknown'] == 'True' else '#AEBBC2' for r in patterns], height=.56)
for i, r in enumerate(patterns):
    left.text(int(r['count']) + 12, i, r['count'], va='center', fontsize=8)
left.set_yticks(range(6), [pattern_short[r['unknown_endpoints']] for r in patterns])
left.invert_yaxis()
left.set_xlim(0, 660)
left.set_xticks([0, 200, 400, 600])
left.set_xlabel('Compounds')
left.spines['left'].set_visible(False)
left.tick_params(axis='y', length=0, pad=4)

mid = fig.add_axes([.50, .23, .19, .58])
m = np.array([float(r['m_min_favouring']) for r in distance if r['attainable'] == 'True'])
bins = np.arange(0, 75, 5)
mid.hist(m, bins=bins, color='#3D6A80', alpha=.9, edgecolor='white', linewidth=.5)
mid.axvline(np.median(m), color='#A17453', linewidth=1.4, linestyle=(0, (4, 3)))
mid.text(2, mid.get_ylim()[1] * .93, f'median {np.median(m):.0f}', color='#A17453', fontsize=8, va='top')
mid.set_xlim(0, 75)
mid.set_xticks([0, 20, 40, 60])
mid.set_yticks([0, 5, 10])
mid.set_xlabel('Exclusive unknowns to resolve')
mid.set_ylabel('Pairs')
mid.grid(axis='y', color='#E4E8EB', linewidth=.5)

right = fig.add_axes([.79, .23, .185, .58])
policies = []
for r in exposure:
    if r['policy'] not in policies:
        policies.append(r['policy'])
for j, f in enumerate(budgets):
    xs, ys = [], []
    for i, p in enumerate(policies):
        row = [r for r in exposure if r['policy'] == p and abs(float(r['budget']) - f) < 1e-9][0]
        xs.append(i + 1)
        ys.append(100 * float(row['unknown_rate']))
    right.scatter(xs, ys, s=16, color=colors[j], label=f'{int(round(f * 100))}%', zorder=3)
right.axhline(27.4, color=muted, linewidth=.9, linestyle=(0, (3, 3)))
right.text(10.4, 26.2, 'pool 27.4%', ha='right', va='top', fontsize=7, color=muted)
right.set_xlim(0, 11)
right.set_ylim(20, 50)
right.set_xticks([1, 5, 10], ['P1', 'P5', 'P10'])
right.set_yticks([20, 30, 40, 50])
right.set_xlabel('Policy')
right.set_ylabel('Unknown rate, %')
right.grid(axis='y', color='#E4E8EB', linewidth=.5)
right.legend(loc='upper left', fontsize=6.5, frameon=False, handletextpad=.2, borderpad=.1)

fig.text(.02, .965, 'a', fontsize=12, weight='bold')
fig.text(.055, .965, 'Unknown endpoints', fontsize=10, weight='bold')
fig.text(.325, .965, 'b', fontsize=12, weight='bold')
fig.text(.35, .965, 'Distance to resolution', fontsize=10, weight='bold')
fig.text(.645, .965, 'c', fontsize=12, weight='bold')
fig.text(.67, .965, 'Selected batches', fontsize=10, weight='bold')
fig.text(.5, .035, 'Panels a and c are descriptive decompositions of the released records, not tests of a missingness\n'
                   'mechanism; resolving exclusive unknowns need not identify a preferred policy.', ha='center', fontsize=8, color=muted, linespacing=1.3)
catalog['figure_4_evidence_gap'] = export(
    fig, 'figure_4_evidence_gap',
    ['joint_unknown_endpoint_patterns.csv', 'identification_distance_primary45.csv', 'selected_unknown_rates.csv'])
catalog['figure_4_evidence_gap']['claim'] = (
    'Unresolved outcomes concentrate in HLM; the distance to identification is 2-70 exclusive unknown records '
    '(median 34); selected batches are enriched in unknown outcomes.')

for stem in ['figure_2_identifiability_map', 'figure_4_evidence_gap']:
    record = catalog[stem]
    record['script'] = 'analysis/plot_v5.py'
    record['interpretation'] = ('Sharp finite-pool identification bounds, not confidence intervals. '
                                'The map and exposure panels are post-review exploratory decompositions of frozen results.')
(figures / 'figure_catalog.json').write_text(json.dumps(catalog, indent=2), encoding='utf-8')
print(json.dumps(catalog, indent=2))
