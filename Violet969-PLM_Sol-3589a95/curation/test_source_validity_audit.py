import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from source_validity_audit import support,supported_sets,bootstrap_indices,binary,multi,read,OUT

class AuditChecks(unittest.TestCase):
 def test_support_deduplicates_and_keeps_clusters_distinct(self):
  d=pd.DataFrame(dict(pfam=['F','F','F','F'],raw_uid=['a','a','b','c'],label=[0,0,1,1],cluster=['X','X','X','Y']))
  r=support(d,['pfam']).iloc[0]
  self.assertEqual((r.n0,r.n1,r.groups0,r.groups1),(1,2,1,2))
  self.assertEqual(supported_sets(d,['pfam'],1),{('F',)})
  self.assertEqual(supported_sets(d,['pfam'],2),set())
 def test_single_class_support(self):
  d=pd.DataFrame(dict(pfam=['F'],raw_uid=['a'],label=[1],cluster=['X']))
  self.assertEqual(support(d,['pfam']).iloc[0].n0,0)
 def test_bootstrap_keeps_group_members(self):
  for ix in bootstrap_indices(['a','a','b','c'],30):
   self.assertEqual(int((ix==0).sum()),int((ix==1).sum()))
 def test_undefined_auc_not_false_chance(self):
  self.assertTrue(np.isnan(binary([0,0],[.1,.2])['auroc']))
  self.assertEqual(binary([0,1],[.5,.5])['auroc'],.5)
 def test_multiclass_absent_validation_classes(self):
  r=multi(np.array([0,1,1]),np.tile([.6,.3,.1],(3,1)),3)
  self.assertEqual(r['represented_classes'],2);self.assertEqual(r['balanced_accuracy'],.5)
 def test_complete_fixed_predictions(self):
  v=read(OUT/'evaluation_manifest.tsv');self.assertEqual(len(v),4000)
  p=read(OUT/'solubility_predictions.tsv')
  for _,g in p.groupby(['protocol','model']):
   self.assertEqual(len(g),4000);self.assertEqual(set(g.raw_uid),set(v.raw_uid))
   self.assertTrue(g.probability.between(0,1).all())
 def test_policy_preserves_all_and_no_inferred_negatives(self):
  p=read(OUT/'policy_manifest.tsv');self.assertEqual(len(p),13880)
  self.assertTrue(p.A_all.all() and p.C_canonical.all() and p.concept_eligible.all())
  self.assertFalse(p.verified_property_eligible.any())
  self.assertFalse(p.archive_explicit_negative_measurement.any())
  self.assertEqual(int(p.archive_explicit_soluble_status.sum()),2)
 def test_recovery_all_samples_accounted_for(self):
  t=read(OUT/'targettrack_sample_matches.tsv');d=read(OUT/'dnasu_sample_matches.tsv')
  self.assertEqual(len(t),24);self.assertTrue(t.recovered.all())
  self.assertEqual(len(d),12);self.assertTrue(d.published_label_agrees.all())
  self.assertFalse(d.exact_listed_protein_match.any())
if __name__=='__main__':unittest.main()
