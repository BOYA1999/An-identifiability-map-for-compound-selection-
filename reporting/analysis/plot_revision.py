import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

revision = Path(__file__).resolve().parent.parent
figures = revision / 'figures'
read = lambda name: list(csv.DictReader((revision / 'tables' / name).open(encoding='utf-8')))
all_pairs = read('paired_decomposition_primary_45.csv')
representatives = read('representative_pairs.csv')
policies = ['rf_paired', 'rf_independent', 'rf_point_min', 'rf_point_sum', 'dmpnn_paired', 'dmpnn_independent', 'dmpnn_point_min', 'dmpnn_point_sum', 'direct_joint_classifier', 'marginal_classifiers_product']
names = ['RF paired', 'RF independent', 'RF point-min', 'RF point-sum', 'DMPNN paired', 'DMPNN independent', 'DMPNN point-min', 'DMPNN point-sum', 'Direct joint classifier', 'Marginal classifier product']
label = dict(zip(policies, names))
codes = dict(zip(policies, [f'P{i + 1}' for i in range(10)]))
colors = ['#3D6A80', '#A17453', '#697D68']
ink, muted = '#253744', '#66727A'
plt.style.use(Path(__file__).parent / 'academic.mplstyle')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.labelsize': 9, 'xtick.labelsize': 8,
                     'ytick.labelsize': 8, 'axes.grid': False, 'axes.edgecolor': '#BBC3C8', 'text.color': ink,
                     'axes.labelcolor': ink, 'figure.dpi': 150, 'savefig.dpi': 300, 'svg.fonttype': 'none'})

def export(fig, stem):
    for extension in ['png', 'pdf', 'svg']:
        options = {'metadata': {'CreationDate': None, 'ModDate': None}} if extension == 'pdf' else {}
        fig.savefig(figures / f'{stem}.{extension}', dpi=300, bbox_inches='tight', pad_inches=.06, **options)
    plt.close(fig)

def interval(ax, row, y, color):
    lo, hi = float(row['lower_pp']), float(row['upper_pp'])
    ax.plot([lo, hi], [y, y], color=color, lw=2.8, solid_capstyle='butt')
    ax.plot([lo, lo], [y - .09, y + .09], color=color, lw=1.2)
    ax.plot([hi, hi], [y - .09, y + .09], color=color, lw=1.2)

fig = plt.figure(figsize=(7.6, 4.8))
ax = fig.add_axes([.35, .47, .60, .43])
for index, row in enumerate(representatives):
    y = 2 - index
    interval(ax, row, y, colors[index])
    ax.text(float(row['lower_pp']) - 1.3, y + .16, f"{float(row['lower_pp']):+.2f}", ha='right', va='center', fontsize=8, color=colors[index])
    ax.text(float(row['upper_pp']) + 1.3, y + .16, f"{float(row['upper_pp']):+.2f}", ha='left', va='center', fontsize=8, color=colors[index])
ax.axvline(0, color='#A9B2B8', lw=.8, linestyle=(0, (3, 3)), zorder=0)
ax.set_xlim(-47, 49)
ax.set_xticks([-40, -20, 0, 20, 40])
ax.set_ylim(-.5, 2.6)
ax.set_yticks([2, 1, 0], [f"{label[r['policy_a']]} −\n{label[r['policy_b']]}" for r in representatives])
ax.tick_params(axis='y', length=0, pad=10)
ax.spines['left'].set_visible(False)
ax.set_xlabel('Joint scenario-attainment difference, A − B (percentage points)', labelpad=7)
fig.text(.02, .97, 'a', weight='bold', fontsize=12)
fig.text(.055, .97, 'Paired identification bounds at the 10% budget (k = 229)', fontsize=10, weight='bold')
fig.text(.02, .345, 'b', weight='bold', fontsize=12)
fig.text(.055, .345, 'What determines each bound', fontsize=10, weight='bold')
tax = fig.add_axes([.02, .07, .96, .24])
tax.set_axis_off()
headers = ['Width rank', 'Shared\n(unknown)', 'Exclusive\nper side', 'Known-success\ndifference, d', 'Unknown\nA-only, uA', 'Unknown\nB-only, uB', 'Width\n(pp)']
cells = [[r['width_rank'] + ' / 45', f"{r['overlap']} ({r['shared_unknown']})", r['exclusive_per_side'],
          f"{int(r['known_success_difference']):+d}", r['unknown_a_only'], r['unknown_b_only'], f"{float(r['width_pp']):.2f}"] for r in representatives]
table = tax.table(cellText=cells, colLabels=headers, cellLoc='center', colLoc='center', bbox=[0, 0, 1, 1], colWidths=[.13, .15, .14, .19, .14, .14, .11])
table.auto_set_font_size(False)
table.set_fontsize(8)
for (r, c), cell in table.get_celld().items():
    cell.set_edgecolor('white')
    cell.set_facecolor('#EDF0F2' if r == 0 else '#F7F8F8' if r % 2 else 'white')
    if r > 0 and c == 0:
        cell.get_text().set_color(colors[r - 1])
        cell.get_text().set_weight('bold')
fig.text(.5, .018, 'Bounds = [(d − uB)/k, (d + uA)/k]; shared candidates cancel. No point estimate is imputed.', ha='center', fontsize=8, color=muted)
export(fig, 'figure_s4_representative_pairs')

fig, axes = plt.subplots(1, 2, figsize=(7.6, 7.7))
fig.subplots_adjust(left=.12, right=.97, bottom=.24, top=.91, wspace=.31)
highlight = {(r['policy_a'], r['policy_b']): colors[i] for i, r in enumerate(representatives)}
for panel, ax in enumerate(axes):
    rows = all_pairs[:23] if panel == 0 else all_pairs[23:]
    for index, row in enumerate(rows):
        interval(ax, row, index, highlight.get((row['policy_a'], row['policy_b']), '#879AA6'))
    ax.set_yticks(range(len(rows)), [f"{codes[r['policy_a']]} − {codes[r['policy_b']]}" for r in rows])
    ax.set_ylim(22.7, -.7)
    ax.set_xlim(-45, 45)
    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.axvline(0, color='#A9B2B8', lw=.8, linestyle=(0, (3, 3)), zorder=0)
    ax.tick_params(axis='y', length=0, pad=5, labelsize=8)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('A − B (percentage points)')
fig.text(.04, .963, 'All 45 later-cohort paired comparisons at the 10% budget', fontsize=11, weight='bold')
fig.text(.04, .934, 'Identification bounds for fixed k = 229 batches; every interval strictly crosses zero.', fontsize=9, color=muted)
for i, (code, name) in enumerate(zip(codes.values(), names)):
    fig.text(.08 + .47 * (i // 5), .176 - .025 * (i % 5), f'{code}  {name}', fontsize=8)
fig.text(.5, .029, 'Saved policy order and orientation are retained. Colors identify the three Figure S4 examples.', ha='center', fontsize=8, color=muted)
export(fig, 'figure_s1_all_primary_contrasts')

catalog_path = figures / 'figure_catalog.json'
catalog = json.loads(catalog_path.read_text()) if catalog_path.exists() else {}
for stem, sources, claim in [
    ('figure_s4_representative_pairs', ['tables/representative_pairs.csv'],
     'Three deterministically chosen pairs illustrate differences in precision; none identifies a strict direction.'),
    ('figure_s1_all_primary_contrasts', ['tables/paired_decomposition_primary_45.csv'],
     'All 45 primary pairs are shown in saved orientation; colors match the Figure S4 representative pairs.')
]:
    catalog[stem] = {'figure': stem, 'surface': 'supplement', 'sources': sources,
        'claim': claim, 'script': 'analysis/plot_revision.py',
        'exports': [f'figures/{stem}.{ext}' for ext in ['png', 'pdf', 'svg']],
        'source_sha256': {name: hashlib.sha256((revision / name).read_bytes()).hexdigest() for name in sources},
        'interpretation': 'Fixed finite-pool identification bounds, not confidence intervals. The selected pairs illustrate different interval widths.'}
catalog_path.write_text(json.dumps(catalog, indent=2), encoding='utf-8')
print('Exported current Figure S1 and Figure S4 as PNG, PDF and SVG.')
