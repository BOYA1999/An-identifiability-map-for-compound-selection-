"""Figure 4: internal positive control, greedy schedule and conditional coverage route."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent.parent
figures = root / 'figures'
read = lambda name: list(csv.DictReader((root / 'tables' / name).open(encoding='utf-8')))
lat = {r['profile']: r for r in read('profile_map_later.csv')}
dev = {r['profile']: r for r in read('profile_map_development.csv')}
frontier = read('measurement_frontier.csv')
cov = read('coverage_requirement_primary45.csv')

order = ['LogD', 'KSOL', 'HLM CLint', 'LogD + KSOL', 'LogD + HLM CLint', 'KSOL + HLM CLint', 'All three']
short = {'LogD': 'LogD', 'KSOL': 'KSOL', 'HLM CLint': 'HLM CLint', 'LogD + KSOL': 'LogD + KSOL',
         'LogD + HLM CLint': 'LogD + HLM', 'KSOL + HLM CLint': 'KSOL + HLM', 'All three': 'All three'}
colors = ['#3D6A80', '#A17453', '#697D68']
ink, muted = '#253744', '#66727A'
plt.style.use(Path(__file__).parent / 'academic.mplstyle')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8, 'axes.grid': False, 'axes.edgecolor': '#BBC3C8', 'text.color': ink,
                     'axes.labelcolor': ink, 'figure.dpi': 150, 'savefig.dpi': 300, 'svg.fonttype': 'none'})

fig = plt.figure(figsize=(7.6, 3.9))
# panel a: development versus later at the 20% budget
axa = fig.add_axes([.13, .22, .25, .60])
y = np.arange(len(order))
dev20 = [int(dev[p]['identified_20pct']) for p in order]
lat20 = [int(lat[p]['identified_20pct']) for p in order]
axa.barh(y - .18, dev20, height=.34, color='#A17453', label='Development (22.6% missing)')
axa.barh(y + .18, lat20, height=.34, color='#879AA6', label='Later (55.7% missing)')
for i, (d, l) in enumerate(zip(dev20, lat20)):
    axa.text(d + .5, i - .18, str(d), va='center', fontsize=7.5, color='#A17453')
    axa.text(l + .5, i + .18, str(l), va='center', fontsize=7.5, color=muted)
axa.set_yticks(y, [short[p] for p in order])
axa.invert_yaxis()
axa.set_xlim(0, 58)
axa.set_xticks([0, 20, 40])
axa.set_xlabel('Identified, of 45')
axa.spines['left'].set_visible(False)
axa.tick_params(axis='y', length=0, pad=4)
axa.set_ylim(6.7, -.6)
axa.scatter([20], [5.55], marker='s', s=22, color='#A17453', clip_on=False)
axa.text(22, 5.55, 'Development (22.6% missing)', fontsize=6.5, va='center')
axa.scatter([20], [6.15], marker='s', s=22, color='#879AA6', clip_on=False)
axa.text(22, 6.15, 'Later (55.7% missing)', fontsize=6.5, va='center')

# panel b: measurement frontier
axb = fig.add_axes([.47, .22, .23, .60])
m = np.array([int(r['measurements']) for r in frontier])
s = np.array([int(r['satisfied']) for r in frontier])
axb.step(m, s, where='post', color='#3D6A80', linewidth=1.6)
marks = {}
for r in frontier:
    marks[int(r['satisfied'])] = int(r['measurements'])
for pairs, offset in [(10, (7, 2)), (30, (7, -10)), (44, (-16, 7))]:
    if pairs in marks:
        axb.scatter([marks[pairs]], [pairs], s=26, color='#A17453', zorder=3, edgecolors='white', linewidth=.5)
        axb.annotate(f'{marks[pairs]}', (marks[pairs], pairs), xytext=offset, textcoords='offset points', fontsize=7.5, color='#A17453')
axb.set_xlim(0, 230)
axb.set_ylim(-2, 48)
axb.set_xticks([0, 50, 100, 150, 200])
axb.set_yticks([0, 10, 20, 30, 40])
axb.set_xlabel('Resolved unknown records')
axb.set_ylabel('Comparisons made identifiable')
axb.grid(axis='y', color='#E4E8EB', linewidth=.5)

# panel c: coverage requirement
axc = fig.add_axes([.79, .22, .19, .60])
cur = [100 * float(r['current_known_fraction']) for r in cov if r['attainable'] == 'True']
req = [100 * float(r['required_known_fraction']) for r in cov if r['attainable'] == 'True']
axc.plot([50, 100], [50, 100], color=muted, linewidth=.8, linestyle=(0, (3, 3)))
axc.scatter(cur, req, s=16, color='#3D6A80', alpha=.8, edgecolors='white', linewidth=.3)
axc.set_xlim(50, 100)
axc.set_ylim(74, 101)
axc.set_xticks([50, 75, 100])
axc.set_yticks([75, 85, 95])
axc.set_xlabel('Known now, %')
axc.set_ylabel('Required on fixed-d route, %')
axc.grid(color='#E4E8EB', linewidth=.4)

fig.text(.02, .955, 'a', fontsize=12, weight='bold')
fig.text(.055, .955, 'Positive control', fontsize=10, weight='bold')
fig.text(.40, .955, 'b', fontsize=12, weight='bold')
fig.text(.428, .955, 'Greedy schedule', fontsize=10, weight='bold')
fig.text(.72, .955, 'c', fontsize=12, weight='bold')
fig.text(.748, .955, 'Fixed-d failure route', fontsize=10, weight='bold')
fig.text(.5, .02, 'Panel a: internal development (clearance missing for 288 of 1,272) and later cohort (1,272 of 2,282).\n'
                  'Panel b: a feasible best-case greedy schedule with all resolved records failing; no proven minimum or expected cost.\n'
                  'Panel c: observed success difference d is fixed and opposite-side unknowns resolve as failures.',
         ha='center', fontsize=7.2, color=muted, linespacing=1.35)

sources = ['profile_map_later.csv', 'profile_map_development.csv', 'measurement_frontier.csv', 'coverage_requirement_primary45.csv']
for extension in ['png', 'pdf', 'svg']:
    options = {'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}
    fig.savefig(figures / f'figure_3_strengthening.{extension}', dpi=300, bbox_inches='tight', pad_inches=.06, **options)
plt.close(fig)
catalog_path = figures / 'figure_catalog.json'
catalog = json.loads(catalog_path.read_text()) if catalog_path.exists() else {}
catalog['figure_3_strengthening'] = {
    'figure': 'figure_3_strengthening', 'sources': sources,
    'exports': [f'figures/figure_3_strengthening.{e}' for e in ['png', 'pdf', 'svg']],
    'source_sha256': {s: hashlib.sha256((root / 'tables' / s).read_bytes()).hexdigest() for s in sources},
    'claim': 'The procedure is not uniformly uninformative in the internal development cohort. A feasible best-case greedy '
             'schedule resolves all 44 nonzero-d later comparisons after 217 measurements, with every resolved record failing. '
             'For all 44 pairs, observed opposite-side coverage is below the threshold for the fixed-d failure route.',
    'script': 'analysis/plot_figure4_v6.py',
    'interpretation': 'The greedy schedule is feasible, not a proven minimum or '
                      'expected cost. The coverage threshold is conditional on fixed observed d and opposite-side failures.'}
catalog_path.write_text(json.dumps(catalog, indent=2), encoding='utf-8')
print('figure 4 exported')
