import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
out = root / 'reporting/figures'
out.mkdir(parents=True, exist_ok=True)
paths = {'development': root / 'artifacts/expansionrx_selection_identification_development_v1', 'later': root / 'artifacts/expansionrx_later_identification_v1'}
tables = {cohort: pd.read_csv(path / 'policy_metrics.csv') for cohort, path in paths.items()}
pairs = {cohort: pd.read_csv(path / 'pairwise_identification.csv') for cohort, path in paths.items()}
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'svg.fonttype': 'none', 'savefig.dpi': 300})
blue, orange, gray, purple = '#0072B2', '#D55E00', '#B8C0C8', '#8B5C9D'
labels = ['RF paired', 'RF independent', 'RF minimum margin', 'RF sum of margins', 'DMPNN paired', 'DMPNN independent', 'DMPNN minimum margin', 'DMPNN sum of margins', 'Joint classifier', 'Product of classifiers']

def finish(fig, name):
    fig.savefig(out / f'{name}.png', bbox_inches='tight', facecolor='white')
    fig.savefig(out / f'{name}.pdf', bbox_inches='tight', facecolor='white')
    plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0), layout='constrained')
development_source = root / 'artifacts/expansionrx_decision_pilot_v1_attempt2/scored_calibration_development.csv'
development = pd.read_csv(development_source).query("subset == 'development'")
later_source = paths['later'] / 'scored_outcomes.csv'
cohorts = [development, pd.read_csv(later_source)]
status = np.array([[int(z.success_lower.sum()), int((z.success_upper == 0).sum()), int((z.success_lower != z.success_upper).sum())] for z in cohorts])
totals = status.sum(axis=1)
left = np.zeros(2)
for j, (label, color) in enumerate(zip(['Success known', 'Failure known', 'Joint outcome unknown'], [blue, gray, orange])):
    widths = status[:, j] / totals
    axes[0].barh([1, 0], widths, left=left, color=color, label=label, height=.5)
    for y, x, width, n in zip([1, 0], left, widths, status[:, j]):
        axes[0].text(x+width/2, y, str(n), ha='center', va='center', fontsize=8, color='white' if j != 1 else 'black')
    left += widths
axes[0].set(yticks=[1, 0], yticklabels=['Development\nN = 1,272', 'Later cohort\nN = 2,282'], xlim=(0, 1), xlabel='Fraction of the full candidate pool', title='a  Decidability of the joint scenario')
axes[0].legend(loc='upper center', bbox_to_anchor=(.5, -.2), frameon=False, fontsize=7, ncol=1)
exact = np.array([int(z.all_exact.sum()) for z in cohorts])
axes[1].barh([1, 0], exact/totals, color=blue, height=.5, label='All three values exact')
axes[1].barh([1, 0], 1-exact/totals, left=exact/totals, color=gray, height=.5, label='At least one value non-exact or absent')
for y, a, n in zip([1, 0], exact, totals):
    axes[1].text(a/n/2, y, str(a), ha='center', va='center', color='white', fontsize=8)
    axes[1].text((1+a/n)/2, y, str(n-a), ha='center', va='center', fontsize=8)
axes[1].set(yticks=[], xlim=(0, 1), xlabel='Fraction of the full candidate pool', title='b  Availability of exact values')
axes[1].legend(loc='upper center', bbox_to_anchor=(.5, -.2), frameon=False, fontsize=7)
finish(fig, 'figure_1_observation_scope')
pd.DataFrame({'cohort': ['development', 'later'], 'known_success': status[:, 0], 'known_failure': status[:, 1], 'unknown': status[:, 2], 'all_exact': exact, 'total': totals}).to_csv(out / 'figure_1_source.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(8.1, 4.8), sharey=True, layout='constrained')
source_rows = []
for ax, (cohort, table), letter in zip(axes, tables.items(), ['a', 'b']):
    z = table.loc[table.fraction.eq(.1)].reset_index(drop=True)
    y = np.arange(10)
    ax.hlines(y, z.success_lower, z.success_upper, color=blue, linewidth=3, label='Full-pool outcome bounds')
    ax.scatter(z.success_lower, y, color=blue, s=12)
    ax.scatter(z.success_upper, y, color=blue, s=12)
    ax.scatter(z.post_exact_success, y-.13, color=orange, marker='^', s=28, label='Exact values within full-pool batch', zorder=3)
    ax.scatter(z.restricted_success, y+.13, color=purple, marker='D', s=19, label='Select after restricting to exact values', zorder=3)
    ax.set(xlim=(0, 1), yticks=y, yticklabels=labels, xlabel='Scenario-success fraction', title=f'{letter}  {cohort.capitalize()}  |  10% budget')
    ax.grid(axis='x', color='#E6E6E6', linewidth=.5)
    source_rows.append(z.assign(cohort=cohort))
axes[0].invert_yaxis()
handles, names = axes[0].get_legend_handles_labels()
fig.legend(handles, names, loc='outside lower center', ncol=1, frameon=False, fontsize=8)
finish(fig, 'figure_s2_policy_success')
pd.concat(source_rows).to_csv(out / 'figure_s2_source.csv', index=False)

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), sharey=True, layout='constrained')
order = ['supported', 'strictly_contradicted', 'unresolved', 'tied', 'unavailable']
colors = [blue, orange, gray, '#4D4D4D', purple]
counts = []
for ax, (cohort, table), letter in zip(axes, pairs.items(), ['a', 'b']):
    x = np.arange(6)
    bottom = np.zeros(6)
    for state, color in zip(order, colors):
        n = np.array([int((table.loc[table.fraction.eq(q), column] == state).sum()) for q in [.05, .1, .2] for column in ['restricted_status', 'post_exact_status']])
        ax.bar(x, n, bottom=bottom, color=color, width=.72, label=state.replace('_', ' ').capitalize())
        for i, val in enumerate(n):
            if val >= 2:
                ax.text(i, bottom[i]+val/2, str(val), ha='center', va='center', fontsize=8, color='black' if state == 'unresolved' else 'white')
        for i, val in enumerate(n):
            counts.append({'cohort': cohort, 'fraction': [.05, .1, .2][i//2], 'report': ['restricted', 'within_batch_exact'][i%2], 'status': state, 'count': int(val)})
        bottom += n
    ax.set(xticks=x, xticklabels=['R\n5%', 'W\n5%', 'R\n10%', 'W\n10%', 'R\n20%', 'W\n20%'], ylim=(0, 49), title=f'{letter}  {cohort.capitalize()}', ylabel='Number of policy pairs (45 total)')
axes[1].set_ylabel('')
handles, names = axes[0].get_legend_handles_labels()
fig.legend(handles[:4], names[:4], loc='outside lower center', ncol=2, frameon=False, fontsize=8)
finish(fig, 'figure_s3_preference_identification')
pd.DataFrame(counts).to_csv(out / 'figure_s3_source.csv', index=False)

sources = [path / name for path in paths.values() for name in ['policy_metrics.csv', 'pairwise_identification.csv']] + [development_source, later_source]
(out / 'source_manifest.json').write_text(json.dumps({str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}, indent=2), encoding='utf-8')
print('Rendered three figures with PDF/PNG exports and numeric source tables.')
