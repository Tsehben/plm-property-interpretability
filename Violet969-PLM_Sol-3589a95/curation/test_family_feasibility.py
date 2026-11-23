"""Synthetic tests of thresholds and family/cluster/source distinctions."""
import unittest
import pandas as pd
from family_feasibility import family_counts, qualify

class FeasibilityTests(unittest.TestCase):
    def fixture(self,n):
        return pd.DataFrame([dict(raw_uid=f'{y}:{i}',pfam='PF00001',label=y,cluster=f'{y}:{i}',center='A' if y==0 else 'B',split='train') for y in (0,1) for i in range(n)])

    def test_threshold_boundaries_and_same_source(self):
        for n in [2,3,4,5]:
            m=self.fixture(n)
            for k in [3,5]:
                self.assertEqual(len(qualify(family_counts(m,['pfam']),k)),int(n>=k))
                self.assertEqual(len(qualify(family_counts(m,['pfam','center']),k)),0)

    def test_proteins_do_not_equal_independent_groups(self):
        m=self.fixture(5);m['cluster']='shared'
        r=family_counts(m,['pfam']).iloc[0]
        self.assertEqual((r.n0,r.n1,r.groups_total,r.groups_shared_labels),(5,5,1,1))
        self.assertEqual((r.groups_exclusive0,r.groups_exclusive1),(0,0))

    def test_multifamily_and_duplicate_memberships(self):
        m=self.fixture(3)
        other=m.copy();other['pfam']='PF00002'
        m=pd.concat([m,other,m.iloc[:1]])
        t=family_counts(m,['pfam'])
        self.assertEqual(len(qualify(t,3)),2)
        self.assertEqual(list(t.n_proteins),[6,6])
        self.assertEqual(m.raw_uid.nunique(),6) # Not twelve independent proteins.

    def test_single_class_families_remain_and_empty_schema(self):
        m=self.fixture(3);m=m[m.label==0]
        t=family_counts(m,['pfam'])
        self.assertEqual(len(t),1)
        self.assertEqual(len(qualify(t,3)),0)
        empty=family_counts(m.iloc[:0],['pfam'])
        self.assertEqual(len(qualify(empty,5)),0)

if __name__=='__main__':unittest.main()
