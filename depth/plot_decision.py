from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

root=Path.cwd()/'analysis_depth/decision'
with (root/'max_regret.csv').open(encoding='utf-8',newline='') as f:
 rows=[r for r in csv.DictReader(f) if r['fraction']=='0.10']
rows.sort(key=lambda r:int(r['max_regret_count']))
names={'rf_paired':'RF: paired residuals','rf_independent':'RF: independent residuals','rf_point_min':'RF: minimum margin','rf_point_sum':'RF: summed margin','dmpnn_paired':'DMPNN: paired residuals','dmpnn_independent':'DMPNN: independent residuals','dmpnn_point_min':'DMPNN: minimum margin','dmpnn_point_sum':'DMPNN: summed margin','direct_joint_classifier':'Direct joint classifier','marginal_classifiers_product':'Product of endpoint classifiers'}
labels=[names[r['policy']] for r in rows]; values=[int(r['max_regret_count']) for r in rows]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(7.1,4.25),layout='constrained')
positions=list(range(len(rows)))
ax.barh(positions,values,height=0.63,color='#567b8c',edgecolor='none')
ax.set_yticks(positions,labels)
ax.invert_yaxis()
ax.set_xlim(0,104)
ax.set_xlabel('Maximum possible shortfall (number of compounds)')
ax.set_title('Worst-case shortfall among the ten candidate strategies\n10% budget (229 compounds); a smaller value does not establish superiority',loc='left',fontsize=9.5,pad=12)
ax.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
ax.xaxis.grid(True,color='#dedede',linewidth=0.6)
ax.set_axisbelow(True)
for p,v in zip(positions,values):ax.text(v+1.2,p,str(v),va='center',ha='left',fontsize=9,color='black')
for spine in ['top','right','left']:ax.spines[spine].set_visible(False)
ax.spines['bottom'].set_color('#777777')
ax.tick_params(axis='y',length=0,pad=6)
ax.tick_params(axis='x',length=3,color='#777777')
fig.savefig(root/'figure_5_decision_shortfall.pdf',bbox_inches='tight')
fig.savefig(root/'figure_5_decision_shortfall.png',dpi=400,bbox_inches='tight')
print([(x['policy'],x['max_regret_count']) for x in rows])

