from pathlib import Path
import zipfile, io, hashlib, json
import pandas as pd
import numpy as np
root=Path.cwd()
P=root/'artifacts/expansionrx_later_identification_v1/scored_outcomes.csv'
O=root/'analysis_depth/mechanism'
O.mkdir(parents=True,exist_ok=True)
def read(path):
 b=path.read_bytes(); return pd.read_csv(io.BytesIO(b)),hashlib.sha256(b).hexdigest()
d,h1=read(P)
ranks,h2=read(root/'artifacts/expansionrx_later_fixed_scores_v2/rankings_and_selection.csv')
pairs,h3=read(root/'reporting/tables/paired_decomposition_all_135.csv')
assert len(d)==2282 and d['Molecule Name'].is_unique
assert len(ranks)==22820 and not ranks.duplicated(['policy','Molecule Name']).any()
d['missing']=d['raw_HLM CLint'].isna() | d['raw_HLM CLint'].astype(str).str.strip().eq('')
d['joint_unknown']=d.success_lower.ne(d.success_upper)
d['stratum']=d.state_0.astype(str)+'/'+d.state_1.astype(str)
assert d['missing'].sum()==1272 and d['joint_unknown'].sum()==625
rows=[]; strata=[]; memberships={}
for frac,col in [(0.05,'top05'),(0.10,'top10'),(0.20,'top20')]:
 union=set()
 for policy,g in ranks.groupby('policy',sort=True):
  ids=set(g.loc[g[col].astype(str).str.lower().eq('true'),'Molecule Name']); union |= ids
  assert len(ids)==int(np.ceil(frac*len(d)))
  memberships[(frac,policy)]=ids
  x=d.assign(selected=d['Molecule Name'].isin(ids)); raw=x.groupby('selected').missing.agg(['size','sum','mean'])
  t=x.groupby(['stratum','selected']).missing.agg(['size','sum','mean']).unstack('selected',fill_value=0)
  common=(t[('size',True)]>0)&(t[('size',False)]>0)
  denom=t.loc[common,('size',True)].sum()+t.loc[common,('size',False)].sum()
  weights=(t.loc[common,('size',True)]+t.loc[common,('size',False)])/denom
  adj_s=(weights*t.loc[common,('mean',True)]).sum(); adj_n=(weights*t.loc[common,('mean',False)]).sum()
  rows.append(dict(fraction=frac,policy=policy,n_selected=len(ids),raw_selected=raw.loc[True,'mean'],raw_not_selected=raw.loc[False,'mean'],common_n=denom,common_selected=t.loc[common,('size',True)].sum(),common_not_selected=t.loc[common,('size',False)].sum(),standardized_selected=adj_s,standardized_not_selected=adj_n,standardized_difference=adj_s-adj_n,common_strata=','.join(t.index[common])))
  for st,row in t.iterrows():
   strata.append(dict(fraction=frac,policy=policy,stratum=st,selected_n=row[('size',True)],selected_missing=int(row[('sum',True)]),not_selected_n=row[('size',False)],not_selected_missing=int(row[('sum',False)]),common_support=bool(common.loc[st])))
 ux=d.assign(selected=d['Molecule Name'].isin(union)); u=ux.groupby('selected').agg(n=('missing','size'),missing=('missing','sum'),joint_unknown=('joint_unknown','sum'))
 print('UNION',frac,u.to_dict('index'))
res=pd.DataFrame(rows); res.to_csv(O/'policy_standardization.csv',index=False)
pd.DataFrame(strata).to_csv(O/'policy_strata.csv',index=False)
q=d.joint_unknown.mean(); out=[]
for _,r in pairs.iterrows():
 k=int(r.k); overlap=int(r.overlap); ex=2*(k-overlap); unk=int(r.exclusive_unknown_total)
 a=memberships[(float(r.fraction),r.policy_a)]; b=memberships[(float(r.fraction),r.policy_b)]
 assert len(a & b)==overlap
 da=d.set_index('Molecule Name').loc[list(a-b)]; db=d.set_index('Molecule Name').loc[list(b-a)]
 assert int(da.joint_unknown.sum())==int(r.unknown_a_only) and int(db.joint_unknown.sum())==int(r.unknown_b_only)
 assert int(da.success_lower.sum()-db.success_lower.sum())==int(r.known_success_difference)
 assert ex==2*int(r.exclusive_per_side)
 assert unk==int(r.unknown_a_only)+int(r.unknown_b_only)
 assert int(r.upper_numerator)-int(r.lower_numerator)==unk
 width=unk/k; assert abs(width*100-r.width_pp)<1e-10
 assert abs((int(r.known_success_difference)-int(r.unknown_b_only))/k*100-r.lower_pp)<1e-10
 assert abs((int(r.known_success_difference)+int(r.unknown_a_only))/k*100-r.upper_pp)<1e-10
 rate=unk/ex if ex else np.nan
 if ex: assert abs(width-2*(1-overlap/k)*q*rate/q)<1e-12
 out.append(dict(fraction=r.fraction,policy_a=r.policy_a,policy_b=r.policy_b,k=k,overlap=overlap,exclusive_count=ex,known_success_difference=int(r.known_success_difference),unknown_exclusive_count=unk,pool_unknown_rate=q,exclusive_unknown_rate=rate,enrichment=rate/q,width_fraction=width,lower_pp=r.lower_pp,upper_pp=r.upper_pp))
pd.DataFrame(out).to_csv(O/'width_decomposition.csv',index=False)
print('POOL',len(d),int(d.missing.sum()),int(d.joint_unknown.sum()))
print('COMPOSITION',d.groupby('stratum').agg(n=('missing','size'),missing=('missing','sum'),joint_unknown=('joint_unknown','sum')).to_string())
for frac,g in res.groupby('fraction'):
 print('STANDARD',frac,'raw diff range', (g.raw_selected-g.raw_not_selected).min(),(g.raw_selected-g.raw_not_selected).max(),'adj diff range',g.standardized_difference.min(),g.standardized_difference.max(),'mean adj',g.standardized_difference.mean(),'common frac',g.common_n.min()/2282,g.common_n.max()/2282)
 print(g[['policy','raw_selected','raw_not_selected','standardized_selected','standardized_not_selected','common_n']].to_string(index=False))
for frac,g in pd.DataFrame(out).groupby('fraction'):
 print('WIDTH',frac,'width median/range',g.width_fraction.median(),g.width_fraction.min(),g.width_fraction.max(),'overlap',g.overlap.min(),g.overlap.max(),'enrichment',g.enrichment.min(),g.enrichment.max())
(O/'input_hashes.json').write_text(json.dumps({'source_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'s8_sha256':h1,'s10_sha256':h2,'s13_sha256':h3},indent=2),encoding='utf8')

union_rows=[]; union_strata=[]
for frac in [0.05,0.1,0.2]:
 ids=set().union(*(v for (f,p),v in memberships.items() if f==frac)); x=d.assign(selected=d['Molecule Name'].isin(ids))
 t=x.groupby(['stratum','selected']).missing.agg(['size','sum','mean']).unstack('selected',fill_value=0)
 common=(t[('size',True)]>0)&(t[('size',False)]>0)
 n=t.loc[common,('size',True)].sum()+t.loc[common,('size',False)].sum()
 w=(t.loc[common,('size',True)]+t.loc[common,('size',False)])/n
 z=x.groupby('selected').agg(n=('missing','size'),missing=('missing','sum'),unknown=('joint_unknown','sum'))
 for selected in [False,True]:
  group=x[x.selected.eq(selected)]; both=group.state_0.eq(1)&group.state_1.eq(1)
  union_rows.append(dict(fraction=frac,selected=selected,n=len(group),hlm_missing=int(group.missing.sum()),joint_unknown=int(group.joint_unknown.sum()),both_pass=int(both.sum()),both_pass_hlm_missing=int(group.loc[both,'missing'].sum()),common_n=n,standardized_hlm_missing=(w*t.loc[common,('mean',selected)]).sum()))
 for s,row in t.iterrows():
  union_strata.append(dict(fraction=frac,stratum=s,selected_n=row[('size',True)],selected_missing=row[('sum',True)],not_selected_n=row[('size',False)],not_selected_missing=row[('sum',False)],common_support=bool(common.loc[s])))
pd.DataFrame(union_rows).to_csv(O/'union_summary.csv',index=False)
pd.DataFrame(union_strata).to_csv(O/'union_strata.csv',index=False)
