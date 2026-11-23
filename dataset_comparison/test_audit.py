import unittest
from audit import *
class AuditTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.r=read(OUT/'records.tsv.gz');cls.m=read(OUT/'policy_groups.tsv.gz');cls.d=read(OUT/'domains.tsv.gz')
 def test_records_and_hashes(self):
  r=self.r;self.assertTrue(r.record_id.is_unique);self.assertTrue((r.sequence.map(sha)==r.sha256).all());self.assertTrue((r.sequence.str.len()==r.length).all())
  self.assertEqual(len(r[r.dataset=='PSI_original']),12216);self.assertEqual(len(r[r.dataset=='PSI_NetSolP']),11226);self.assertEqual(len(r[r.dataset=='NESG_raw']),9703)
 def test_frozen_uesolds(self):
  r=self.r[self.r.dataset=='UESolDS_development'];original=read(CUR/'results/cohort_manifest.tsv');o=original[original.clean&original.cached].set_index('raw_uid');self.assertEqual(len(r),13880)
  for col in ['sha256','label','split','cluster','center']:
   a=r.set_index('raw_uid')[col].sort_index();b=o[col].sort_index()
   if col=='label':a=pd.to_numeric(a);b=pd.to_numeric(b)
   self.assertTrue(a.index.equals(b.index));self.assertTrue((a==b).all(),col)
 def test_previous_inputs_preserved(self):
  for p,h in json.loads((OUT/'input_checksums.json').read_text()).items():self.assertEqual(digest(ROOT.parent/p),h,p)
 def test_reference_pins(self):
  for ref in json.loads((ROOT/'reference_manifest.json').read_text()):self.assertEqual(digest(ref['local_path']),ref['local_sha256'])
 def test_score_mapping(self):
  d=read(OUT/'NetSolP_NESG_score_checks.tsv');self.assertEqual(len(d),1323);self.assertTrue(d[['exact_sequence','E_agrees','S_agrees']].all().all());self.assertEqual(set(d.loc[~d.label_rule_agrees,'sid']),{'MAR120','MAR123','PFR1'})
  raw=self.r.set_index('record_id');q=self.m[self.m.policy=='NESG_unique_E45_S45_vs0'];a=raw.loc[q.record_id];self.assertTrue(pd.to_numeric(a.exp).isin([4,5]).all());self.assertTrue(pd.to_numeric(a.sol).isin([0,4,5]).all());self.assertEqual(len(q),3762)
 def test_pfam_thresholds_and_coordinates(self):
  d=self.d;self.assertTrue((d.sequence_score>=d.sequence_ga).all());self.assertTrue((d.domain_score>=d.domain_ga).all());self.assertTrue(((d.ali_from>=1)&(d.ali_to<=d.sequence_length)&(d.ali_from<=d.ali_to)).all())
  m=read(OUT/'family_memberships.tsv.gz');self.assertFalse(m.duplicated(['sha256','pfam']).any())
 def test_clusters_and_partitions(self):
  groups=read(OUT/'sequence_groups.tsv.gz').set_index('sha256');e=pd.read_csv(OUT/'homology_edges.tsv.gz',sep='\t',header=None,usecols=[0,1]);self.assertTrue((e[0].map(groups.audit_cluster)==e[1].map(groups.audit_cluster)).all());self.assertEqual(groups.groupby('audit_cluster').audit_fold.nunique().max(),1)
 def test_no_double_count_representatives(self):
  d=read(OUT/'representative_distributions.tsv');g=d.groupby(['policy','level','label']);self.assertTrue(np.allclose(g.percent.sum(),100));self.assertTrue((g.N.sum()==g.denominator.first()).all())
 def test_psi_construct_links(self):
  d=read(OUT/'PSI_construct_mapping.tsv');self.assertTrue(d.label_agrees.all());self.assertTrue(d.contained.all());self.assertEqual(set(zip(d.N_removed,d.C_removed)),{(10,0),(0,8)})
if __name__=='__main__':unittest.main()
