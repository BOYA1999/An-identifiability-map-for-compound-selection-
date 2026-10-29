from pathlib import Path
import zipfile, csv, io, json, hashlib, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import csc_matrix

root=Path.cwd()
zpath=root/'artifacts/expansionrx_later_identification_v1/scored_outcomes.csv'
out=root/'analysis_depth/decision'
out.mkdir(exist_ok=True,parents=True)
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
data=read(zpath)
ranks=read(root/'artifacts/expansionrx_later_fixed_scores_v2/rankings_and_selection.csv')
pairs=read(root/'reporting/tables/paired_decomposition_all_135.csv')
print('loaded',len(data),len(ranks),len(pairs),flush=True)
ids=[r['Molecule Name'] for r in data]; assert len(ids)==len(set(ids))==2282
lo={r['Molecule Name']:int(float(r['success_lower'])) for r in data}; hi={r['Molecule Name']:int(float(r['success_upper'])) for r in data}
assert sum(lo[x]!=hi[x] for x in ids)==625
policies=list(dict.fromkeys(r['policy'] for r in ranks)); assert len(policies)==10
allresults=[]; regrets=[]; witnesses={}; audit=[]
for frac,col,k in [('0.05','top05',115),('0.10','top10',229),('0.20','top20',457)]:
    sel={p:{r['Molecule Name'] for r in ranks if r['policy']==p and r[col]=='True'} for p in policies}
    assert all(len(v)==k for v in sel.values())
    known={p:sum(lo[x] for x in sel[p]) for p in policies}
    unknown=[x for x in ids if lo[x]!=hi[x]]
    groups={}
    for x in unknown:
        signature=tuple(int(x in sel[p]) for p in policies)
        if 0<sum(signature)<10: groups.setdefault(signature,[]).append(x)
    signatures=list(groups); mat=np.array(signatures,dtype=int).T
    ns=np.array([len(groups[s]) for s in signatures],dtype=int)
    print('budget',frac,'groups',len(ns),flush=True)
    saved={(r['policy_a'],r['policy_b']):r for r in pairs if float(r['fraction'])==float(frac)}
    assert len(saved)==45
    upper={}
    for a,b in itertools.combinations(policies,2):
        A,B=sel[a]-sel[b],sel[b]-sel[a]
        L=sum(lo[x] for x in A)-sum(hi[x] for x in B)
        U=sum(hi[x] for x in A)-sum(lo[x] for x in B)
        row=saved.get((a,b)); rev=False
        if row is None: row=saved[(b,a)]; rev=True
        expected=(int(row['lower_numerator']),int(row['upper_numerator']))
        if rev: expected=(-expected[1],-expected[0])
        assert (L,U)==expected,(frac,a,b,L,U,expected)
        upper[a,b]=U; upper[b,a]=-L
        audit.append({'fraction':frac,'a':a,'b':b,'L':L,'U':U,'saved_match':True})
    for a in policies:
        rival=max(policies,key=lambda b: upper.get((b,a),0))
        regret=upper.get((rival,a),0)
        regrets.append({'fraction':frac,'k':k,'policy':a,'max_regret_count':regret,'max_regret_percentage_points':100*regret/k,'worst_rival':rival})
        ai=policies.index(a); competitors=[p for p in policies if p!=a]
        D=np.array([mat[ai]-mat[policies.index(b)] for b in competitors])
        constants=np.array([known[a]-known[b] for b in competitors])
        for strict in [0,1]:
            con=LinearConstraint(csc_matrix(D),lb=np.full(9,strict)-constants,ub=np.full(9,np.inf))
            res=milp(c=np.zeros(len(ns)),integrality=np.ones(len(ns)),bounds=Bounds(np.zeros(len(ns)),ns),constraints=con,options={'time_limit':60})
            feasible=res.status==0
            if feasible:
                v=np.rint(res.x).astype(int); assert np.all(v>=0) and np.all(v<=ns)
                scores=np.array([known[p] for p in policies])+mat@v
                margin=min(scores[ai]-scores[policies.index(b)] for b in competitors)
                assert margin>=strict
                witness=[{'signature':''.join(map(str,s)),'success_count':int(v[i]),'members':groups[s]} for i,s in enumerate(signatures) if v[i]]
                witnesses[f'{frac}|{a}|{strict}']={'scores':dict(zip(policies,map(int,scores))),'groups':witness}
            else: margin=None
            allresults.append({'fraction':frac,'k':k,'policy':a,'strict':strict,'status':res.status,'solver_message':res.message,'feasible':feasible if res.status in (0,2) else None,'verified_margin':margin})
            print(frac,a,'strict',strict,'status',res.status,'margin',margin,flush=True)
def write(name,rows):
    with (out/name).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader();w.writerows(rows)
write('possible_optimal.csv',allresults);write('max_regret.csv',regrets);write('pair_check.csv',audit)
(out/'witnesses.json').write_text(json.dumps(witnesses,ensure_ascii=False,indent=2),encoding='utf-8')
info={'source_file':zpath.relative_to(root).as_posix(),'source_sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'policies':policies,'note':'Each signature count specifies first success_count input-order IDs as successes; all other unknowns 0. This creates one common outcome assignment per witness; excluded all/none signatures fixed 0. Optimization only uses existing finite pool and outcomes. Solver status 2=infeasible; other nonzero not treated as infeasible.'}
(out/'metadata.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
print('regrets',json.dumps(regrets),flush=True)
