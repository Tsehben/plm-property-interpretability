"""Evaluation-only tests; synthetic prediction fixtures never become research results."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from property_inputs import load_inputs,CONTINUOUS
from within_property_evaluation import make_strata,plan,read_predictions,evaluate_layer,render
from run_within_property import select_source,predictions_from_archive
from run_layerwise_auroc import make_bundle
from build_within_property_notebook import build


def fixture():
    n=60
    d=pd.DataFrame(dict(record_id=[f"fixture:{i}" for i in range(n)],sha256=[f"hash{i}" for i in range(n)],
                        split=["train"]*40+["validation"]*20,label=np.arange(n)%2,audit_cluster=np.arange(n)))
    p=pd.DataFrame({name:np.linspace(-1,1,n)[np.random.default_rng(i).permutation(n)] for i,name in enumerate(CONTINUOUS)})
    m=pd.DataFrame(dict(record_id=d.record_id,pfam="PF_fixture",name="Fixture",competing_overlap=np.arange(n)%3==0))
    c=pd.DataFrame(dict(pfam=["PF_fixture"],eligible=[False]))
    return d,p,m,c


class WithinTests(unittest.TestCase):
    def test_training_edges_do_not_use_validation(self):
        d,p,m,c=load_inputs()
        before=make_strata(d,p,m)[1]
        changed=p.copy(); changed.loc[d.split.eq("validation"),CONTINUOUS]=1e9
        after=make_strata(d,changed,m)[1]
        self.assertEqual(before,after)
        self.assertEqual(set(d.split),{"train","validation"})

    def test_partition_boundaries_missingness_and_multidomain_counting(self):
        d,p,m,c=fixture()
        p.loc[45,"hydrophobicity"]=np.nan
        m=pd.concat([m,m.iloc[:3]],ignore_index=True)
        strata,edges,f=make_strata(d,p,m)
        for control in CONTINUOUS:
            masks=[s["mask"] for s in strata if s["kind"]=="range" and s["control"]==control]
            np.testing.assert_array_equal(np.sum(masks,axis=0),p[control].notna().astype(int))
        self.assertEqual(int(f.loc[f.variant.eq("all_ga"),"validation_n"].iloc[0]),20)
        self.assertLess(int(f.loc[f.variant.eq("noncompeting"),"validation_n"].iloc[0]),20)
        # Equal-valued descriptors form one range, not artificial empty quintiles.
        p["hydrophobicity"]=3
        self.assertEqual(make_strata(d,p,m)[1]["hydrophobicity"],[])

    def test_predictions_align_by_id_and_reject_corruption(self):
        d,p,m,c=fixture(); y=d.label.to_numpy(float)
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/"p.tsv"
            table=d.loc[d.split.eq("validation"),["record_id","sha256","split"]].copy()
            table["prediction"]=np.arange(20)
            table.sample(frac=1,random_state=2).to_csv(path,sep="\t",index=False)
            result=read_predictions(path,d,"solubility",y)
            np.testing.assert_array_equal(result[40:],np.arange(20))
            table.loc[40,"sha256"]="wrong"
            table.to_csv(path,sep="\t",index=False)
            with self.assertRaises(AssertionError): read_predictions(path,d,"solubility",y)

    def test_metrics_use_same_subset_and_no_fitting(self):
        d,p,m,c=fixture()
        targets,jobs,support,edges,families=plan(d,p,m,c)
        rng=np.random.default_rng(4)
        scores={name:rng.normal(size=len(d)) for name,_,_ in targets}
        baseline={name:rng.normal(size=len(d)) for name,_,_ in targets}
        weights=np.ones((20,20),int)
        with patch("sklearn.linear_model.LogisticRegression.fit",side_effect=AssertionError("No fitting allowed")),patch("sklearn.linear_model.Ridge.fit",side_effect=AssertionError("No fitting allowed")):
            result=evaluate_layer(1,scores,baseline,jobs,weights)
        row=result[(result.key=="overall") & (result.target=="solubility") & (result.metric=="auroc")].iloc[0]
        expected=roc_auc_score(d.label.to_numpy()[40:],scores["solubility"][40:])
        comp=roc_auc_score(d.label.to_numpy()[40:],baseline["solubility"][40:])
        self.assertAlmostEqual(row.value,expected)
        self.assertAlmostEqual(row.delta_vs_composition,expected-comp)
        family_jobs=[j for j in jobs if j["info"]["kind"]=="family"]
        self.assertTrue(all(j["info"]["target"]!="PF_fixture" for j in family_jobs))

    def test_incomplete_source_stops_without_work(self):
        class Store:
            names={"esm_stage1/uncontrolled_property_probing_v1/fixture/layer_01/completed.json"}
            def get(self,*args): raise AssertionError("No read needed")
        with self.assertRaisesRegex(RuntimeError,"finish through layer 33"):
            select_source(Store())

    def test_archive_extracts_only_predictions(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); src=root/"src"; (src/"solubility").mkdir(parents=True)
            (src/"solubility/validation_predictions.tsv").write_text("record_id\tprediction\nfixture\t0.2\n")
            np.savez(src/"solubility/parameters.npz",coef=np.ones(2))
            raw=make_bundle(src); digest=hashlib.sha256(raw).hexdigest()
            predictions_from_archive(raw,digest,["solubility/validation_predictions.tsv"],root/"out")
            self.assertTrue((root/"out/solubility/validation_predictions.tsv").exists())
            self.assertFalse((root/"out/solubility/parameters.npz").exists())
            with self.assertRaises(AssertionError): predictions_from_archive(raw,"wrong",[],root/"bad")

    def test_notebook_compiles_and_has_no_extraction_or_fit_commands(self):
        import nbformat
        nb=build(); nbformat.validate(nb)
        for cell in nb.cells:
            if cell.cell_type!="code": continue
            compile(cell.source,"<cell>","exec")
            for node in ast.walk(ast.parse(cell.source)):
                if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=="run_stage1":
                    if node.args and isinstance(node.args[0],ast.Constant) and node.args[0].value=="-c":
                        compile(ast.literal_eval(node.args[1]),"<inner-script>","exec")
        source="\n".join(c.source for c in nb.cells if c.cell_type=="code")
        self.assertNotIn("--activations",source)
        self.assertNotIn("run_stage1('run_property_probes.py'",source)
        self.assertIn("run_within_property.py",source)

    def test_end_to_end_render_on_temporary_fixtures(self):
        d,p,m,c=fixture()
        targets,jobs,support,edges,families=plan(d,p,m,c)
        scores={name:np.linspace(-1,1,len(d)) for name,_,_ in targets}
        metrics=evaluate_layer(0,scores,scores,jobs,np.ones((20,20),int))
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            for layer in range(34):
                metrics.assign(layer=layer).to_csv(folder/f"layer_{layer:02d}.tsv",sep="\t",index=False)
            render(folder,support,c)
            self.assertEqual(len(list((folder/"figures").glob("*.png"))),8)
            table=pd.read_csv(folder/"within_property_metrics.tsv",sep="\t")
            self.assertTrue(table.delta_vs_composition.dropna().eq(0).all())
            self.assertEqual(set(table.layer),set(range(34)))


if __name__=="__main__": unittest.main()
