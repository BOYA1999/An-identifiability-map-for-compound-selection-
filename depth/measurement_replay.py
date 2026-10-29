from pathlib import Path
from itertools import combinations
import hashlib, json, math
import numpy as np
import pandas as pd

root = Path.cwd()
out = root / 'analysis_depth/replay'
out.mkdir(parents=True,exist_ok=True)
paths = [root / 'artifacts/expansionrx_decision_support_v1/adjudicated_training.csv', root / 'artifacts/expansionrx_decision_pilot_v1_attempt2/scored_calibration_development.csv']
a, s = (pd.read_csv(p) for p in paths)
a = a[a.subset == 'development'].copy()
s = s[s.subset == 'development'].copy()
assert len(a) == len(s) == 1272
assert a['Molecule Name'].is_unique and s['Molecule Name'].is_unique
assert set(a['Molecule Name']) == set(s['Molecule Name'])
a = a[a.all_exact].copy()
assert len(a) == 928
s = s.drop(columns=['subset', 'all_exact', 'split_group', 'success_lower', 'success_upper'])
d = a.merge(s, on='Molecule Name', validate='one_to_one').sort_values('Molecule Name').reset_index(drop=True)
assert len(d) == 928
states = d[['state_0','state_1','state_2']].to_numpy(dtype=int)
assert np.isin(states,[0,1]).all()
truth = (states == 1).all(axis=1).astype(int)
assert np.array_equal(truth,d.success_lower.to_numpy()) and np.array_equal(truth,d.success_upper.to_numpy())
raw = d[['LogD','KSOL','HLM CLint']].apply(pd.to_numeric, errors='raise').to_numpy()
assert np.array_equal(states[:,0],((raw[:,0]>=1)&(raw[:,0]<=3)).astype(int))
assert np.array_equal(states[:,1],(raw[:,1]>50).astype(int))
assert np.array_equal(states[:,2],(raw[:,2]<10).astype(int))
rules=['rf_paired','rf_independent','rf_point_min','rf_point_sum','dmpnn_paired','dmpnn_independent','dmpnn_point_min','dmpnn_point_sum','direct_joint_classifier','marginal_classifiers_product']
n=len(d); k=math.ceil(.1*n)
assert k==93 and np.isfinite(d[rules].to_numpy(dtype=float)).all()
sel=np.zeros((10,n), dtype=int)
for i,rule in enumerate(rules):
    order=np.lexsort((d['Molecule Name'].to_numpy(),-d[rule].to_numpy()))
    sel[i,order[:k]]=1
assert (sel.sum(axis=1)==k).all()
pairs=list(combinations(range(10),2))
diffs=np.array([sel[i]-sel[j] for i,j in pairs])
degree=(diffs!=0).sum(axis=0)
actual=diffs@truth
observed_successes=sel@truth
rng=np.random.default_rng(43187)
rows=[]; masks=[]; tests=0
for rep in range(200):
    masked=rng.choice(n, size=math.floor(.4*n), replace=False)
    candidates=masked[(states[masked,0]==1)&(states[masked,1]==1)]
    random_order=rng.permutation(candidates)
    focused_order=candidates[np.lexsort((d['Molecule Name'].to_numpy()[candidates],-degree[candidates]))]
    maskflag=np.zeros(n,dtype=bool); maskflag[candidates]=True
    for index in candidates:
        masks.append({'replicate':rep,'molecule':d.iloc[index]['Molecule Name'],'degree':int(degree[index]),'outcome':int(truth[index]),'random_rank':int(np.flatnonzero(random_order==index)[0]+1),'degree_rank':int(np.flatnonzero(focused_order==index)[0]+1)})
    relevant_order=np.concatenate((random_order[degree[random_order]>0], random_order[degree[random_order]==0]))
    for strategy,order in [('random',random_order),('relevant_random',relevant_order),('degree',focused_order)]:
        for budget in [0,10,25,50,100]:
            reveal=order[:budget]
            unknown=maskflag.copy(); unknown[reveal]=False
            low=truth.copy(); high=truth.copy(); low[unknown]=0; high[unknown]=1
            lower=np.maximum(diffs,0)@low+np.minimum(diffs,0)@high
            upper=np.maximum(diffs,0)@high+np.minimum(diffs,0)@low
            assert np.all(lower<=actual) and np.all(actual<=upper)
            assert np.array_equal(upper-lower,np.abs(diffs)@unknown.astype(int))
            orient=np.zeros((10,10),dtype=int)
            for z,(i,j) in enumerate(pairs):
                orient[i,j]=upper[z]; orient[j,i]=-lower[z]
            regret=orient.max(axis=0)
            realized=observed_successes.max()-observed_successes
            assert np.all(regret>=realized)
            identified=((lower>0)|(upper<0))
            assert np.all((actual[identified]>0)==(lower[identified]>0))
            tests+=1
            rows.append({'replicate':rep,'strategy':strategy,'budget_requested':budget,'budget_actual':len(reveal),'candidate_count':len(candidates),'identified':int(identified.sum()),'width_sum_counts':int((upper-lower).sum()),'minimax_regret_count':int(regret.min()),'revealed_successes':int(truth[reveal].sum())})
res=pd.DataFrame(rows)
res.to_csv(out/'replicates.csv',index=False)
pd.DataFrame(masks).to_csv(out/'mask_and_order.csv',index=False)
summary=(res.groupby(['strategy','budget_requested'],as_index=False).agg(n=('replicate','size'),candidate_min=('candidate_count','min'),candidate_max=('candidate_count','max'),budget_min=('budget_actual','min'),identified_mean=('identified','mean'),identified_median=('identified','median'),width_sum_mean=('width_sum_counts','mean'),regret_mean=('minimax_regret_count','mean'),reveal_success_mean=('revealed_successes','mean')))
summary.to_csv(out/'summary.csv',index=False)
pivot=res.pivot(index=['replicate','budget_requested'],columns='strategy',values=['identified','width_sum_counts','minimax_regret_count'])
paired=pd.DataFrame({metric:pivot[metric]['degree']-pivot[metric]['random'] for metric in ['identified','width_sum_counts','minimax_regret_count']}).reset_index()
paired.to_csv(out/'paired_difference.csv',index=False)
print(summary.to_string(index=False))
print('paired means degree minus random')
print(paired.groupby('budget_requested').mean().to_string())
manifest={'inputs':{p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'n_development':1272,'n_complete_exact':n,'k':k,'seed':43187,'replicates':200,'hidden_HLM_count_per_run':math.floor(.4*n),'actual_successes_per_policy':dict(zip(rules,map(int,observed_successes))),'validation_evaluations':tests,'candidate_min':int(res.candidate_count.min()),'candidate_max':int(res.candidate_count.max())}
(out/'checks.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')

