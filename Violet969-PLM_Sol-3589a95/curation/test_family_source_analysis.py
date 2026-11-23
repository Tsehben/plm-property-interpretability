import unittest
import pandas as pd
from family_source_analysis import decompose,count_groups

class VariationTests(unittest.TestCase):
    def make(self,kind):
        rows=[]
        for f in ['F0','F1']:
            for c in ['C0','C1']:
                for i in range(2):
                    y=int(f=='F1') if kind=='family' else int(c=='C1') if kind=='source' else i
                    rows.append(dict(raw_uid=f+c+str(i),pfam=f,center=c,label=y,cluster=f+c+str(i)))
        return pd.DataFrame(rows)
    def test_between_family(self):
        r=decompose(self.make('family'),'pfam','center')
        self.assertAlmostEqual(r['between_first_fraction'],1)
        self.assertAlmostEqual(r['within_both_fraction'],0)
    def test_source_within_family(self):
        r=decompose(self.make('source'),'pfam','center')
        self.assertAlmostEqual(r['within_first_between_second_fraction'],1)
    def test_within_both(self):
        r=decompose(self.make('within'),'pfam','center')
        self.assertAlmostEqual(r['within_both_fraction'],1)
    def test_multifamily_weighting_and_groups(self):
        m=self.make('within');extra=m.iloc[:1].copy();extra['pfam']='F2'
        m=pd.concat([m,extra]);r=decompose(m,'pfam','center')
        self.assertEqual(r['n_unique'],8)
        self.assertAlmostEqual(r['positive_prevalence'],.5)
        m['cluster']='ONE';t=count_groups(m,['pfam'])
        self.assertEqual(t.iloc[0].groups0,1)
        self.assertEqual(t.iloc[0].groups1,1)
        self.assertEqual(t.iloc[0].exclusive_groups0,0)
class SavedResultTests(unittest.TestCase):
    def test_saved_counts_predictions_and_input_integrity(self):
        import json,hashlib
        import numpy as np
        from pathlib import Path
        from sklearn.metrics import roc_auc_score
        p=Path(__file__).resolve().parent/'family_source_results'
        def read(name):return pd.read_csv(p/name,sep='\t',float_precision='round_trip')
        c=read('focus_69_family_source_counts.tsv');f=read('focus_69_family_counts.tsv')
        a=c.groupby(['split','pfam'])[['n0','n1']].sum().sort_index()
        b=f.set_index(['split','pfam'])[['n0','n1']].sort_index()
        np.testing.assert_array_equal(a.values,b.values)
        r=read('baseline_group_results.tsv');pred=read('baseline_predictions.tsv')
        for (protocol,policy,model),g in pred.groupby(['protocol','policy','model']):
            self.assertTrue(g.raw_uid.is_unique)
            target=r[(r.protocol==protocol)&(r.policy==policy)&(r.model==model)&(r.group_type=='pooled')].iloc[0]
            self.assertAlmostEqual(roc_auc_score(g.label,g.prediction),target.auroc,places=12)
            if model=='source':self.assertTrue((g.groupby('center').prediction.nunique()==1).all())
        rr=r[(r.protocol=='published_split')&(r.policy=='GA_all')&(r.group_type=='exact_family_set_source')]
        np.testing.assert_allclose(rr.auroc,.5)
        meta=json.loads((p/'provenance.json').read_text())
        for name,h in meta['inputs'].items():self.assertEqual(hashlib.sha256((p.parent/name).read_bytes()).hexdigest(),h)

if __name__=='__main__':unittest.main()
