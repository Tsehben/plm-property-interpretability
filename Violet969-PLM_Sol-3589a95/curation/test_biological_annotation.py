"""Scientific edge cases and published external fixtures; unittest, no network."""
import unittest
import hashlib,json
from pathlib import Path
import numpy as np
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from biological_features import *

class Features(unittest.TestCase):
    def test_published_aggrescan_sh3_example(self):
        # 2007 Additional file 3, wild-type SH3. Independent published fixture.
        s='MDETGKELVLALYDYQEKSPREVTMKKGDILTLLNSTNKDWWKVEVNDRQGFVPAAYVKKLD'
        self.assertEqual(len(s),62)
        p,rr,status=aggregation(s)
        self.assertEqual(rr,[(7,14),(30,35),(55,59)])
        self.assertAlmostEqual(p.mean()*100,-13.9,places=1)
        np.testing.assert_allclose(p[2:6],[-.606,-.975,-.890,-.331],atol=.001)
        # The supplementary output disagrees with the prose terminal ordering.
        pe,_,_=aggregation(s,terminal_policy='published_example')
        np.testing.assert_allclose(pe[:6],[-.824,-.824,-.606,-.975,-.890,-.331],atol=.001)
        self.assertAlmostEqual(pe[-1],-.681,places=3)
        self.assertGreater(abs(p[0]-pe[0]),.1)
        mutant=s[:46]+'A'+s[47:]
        p2,rr2,_=aggregation(mutant)
        self.assertEqual(len(rr2),4)
        self.assertAlmostEqual(p2.mean()*100,-11.9,places=1)

    def test_unknowns_not_neutral_or_imputed(self):
        b=biochemistry('AKXDE')
        self.assertTrue(np.isnan(b['net_charge_ph7p0']))
        self.assertTrue(np.isnan(b['hydropathy_mean']))
        self.assertAlmostEqual(sum(b['frac_'+a] for a in AA),.8)
        p,rr,st=aggregation('AKXDE');self.assertTrue(np.isnan(p).all());self.assertEqual(st,'unsupported_input')

    def test_short_windows_not_zero(self):
        self.assertEqual(len(window_values([1.,2.],9)),0)
        b=biochemistry('AK');self.assertTrue(np.isnan(b['hydropathy_w9_max']))
        self.assertEqual(aggregation('AA')[2],'unsupported_input')

    def test_charge_actual_termini_and_pi(self):
        s='MACKRHDEYQW'
        for ph in [7.,7.4]:self.assertAlmostEqual(charge_profile(s,ph).sum(),IsoelectricPoint(s).charge_at_pH(ph),places=12)
        self.assertNotAlmostEqual(charge_profile(s,7.)[:5].sum(),charge_profile(s[:5],7.).sum())
        b=biochemistry('D'*50);self.assertLess(b['pi_estimate'],4.05)
        self.assertLess(abs(IsoelectricPoint('D'*50).charge_at_pH(b['pi_estimate'])),.001)

    def test_coordinates_and_visible_region_qualification(self):
        self.assertEqual(regions([True,True,False,True]),[(1,2),(4,4)])
        self.assertEqual(clip_regions([(1021,1030),(1040,1050)],1022),[(1021,1022)])
        # A full-length qualifying APR remains a partial visible region of 2 aa.
        self.assertEqual(regional_summary(np.ones(1022),[(1021,1022)])['n_regions'],1)
        for n in [5,75,76,175,176,300,301,2001]:
            p,rr,_=aggregation('V'*n);self.assertEqual(len(p),n)
            self.assertTrue(all(1<=a<=b<=n for a,b in rr))

    def test_proline_breaks_aprs_and_entropy(self):
        _,rr,_=aggregation('V'*10+'P'+'V'*10)
        self.assertTrue(all(not a<=11<=b for a,b in rr))
        self.assertAlmostEqual(biochemistry(AA)['entropy_bits'],np.log2(20))
        self.assertAlmostEqual(biochemistry('AAAA')['entropy_bits'],0.)

    def test_tm_regions_are_distinct_from_membrane_association(self):
        seg=[('signal',1,10),('TMhelix',20,39),('reentrant',70,79),('Beta sheet',90,99)]
        a=membrane_summary(seg,100)
        self.assertEqual(a['tm_segment_count'],2)
        self.assertAlmostEqual(a['tm_fraction'],.3)
        self.assertAlmostEqual(a['membrane_associated_fraction'],.4)
        self.assertEqual(membrane_summary(seg,25)['tm_segment_count'],1)
        self.assertAlmostEqual(membrane_summary(seg,25)['tm_fraction'],6/25)

    def test_export_packets_preserve_every_sequence(self):
        from annotate_biology import OUT,load,read
        manifest=OUT/'server_upload_manifest.tsv'
        if not manifest.exists():self.skipTest('Export packets not generated yet')
        d,seqs=load();m=read(manifest)
        nsp=m[m.tool=='NetSurfP3.0']
        self.assertEqual(set(nsp.raw_uid),set(d.raw_uid))
        self.assertTrue(nsp.raw_uid.is_unique)
        signal=m[m.tool=='SignalP6.0'];self.assertEqual(set(signal.raw_uid),set(d.raw_uid));self.assertTrue(signal.raw_uid.is_unique)
        export=json.loads((OUT/'server_upload_export.json').read_text());folder=Path(export['destination'])
        if not folder.exists():self.skipTest('Temporary FASTA packet not present')
        lookup=d.set_index('raw_uid')
        for name,g in m.groupby('file'):
            lines=(folder/name).read_text().splitlines();self.assertEqual(len(lines),2*len(g))
            for r,h,s in zip(g.itertuples(),lines[::2],lines[1::2]):
                self.assertEqual(h,'>'+r.upload_id);self.assertEqual(s,seqs[r.raw_uid])
                self.assertEqual(hashlib.sha256(s.encode()).hexdigest(),lookup.loc[r.raw_uid,'sha256'])

if __name__=='__main__':unittest.main()
