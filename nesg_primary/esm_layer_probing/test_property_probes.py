"""Numerical/leakage/resume tests; synthetic fixtures are temporary, never research outputs."""
import ast
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from property_probes import (RidgeBasis, regression_metrics, binary_metrics, target_list,
                             fit_parameters, fit_target, valid_checkpoint, reuse_solubility)
from property_inputs import load_inputs, family_catalog, SPEC
from run_layerwise_auroc import make_bundle, unpack_bundle
from run_property_probes import local_complete, mark_local, PrivateStore, compatible_solubility
from build_property_notebook import build
from common import HERE, composition


class PropertyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data, cls.props, cls.members, cls.catalog = load_inputs()

    def test_frozen_inputs_and_missingness(self):
        d,p,m,c = self.data,self.props,self.members,self.catalog
        self.assertEqual(d.split.value_counts().to_dict(),{"train":2361,"validation":691})
        self.assertEqual(int(p.apr_fraction.isna().sum()),3)
        self.assertGreater(int(p.apr_fraction.eq(0).sum()),0)
        self.assertTrue(p.log_length.notna().all())
        self.assertEqual(int(c.eligible.sum()),21)
        self.assertTrue(set(m.record_id) <= set(d.record_id))
        self.assertEqual(len(target_list(d,p,m,c)),27)

    def test_training_family_selection_ignores_validation(self):
        # Delete all validation memberships and invert validation solubility labels.
        d=self.data.copy()
        d.loc[d.split.eq("validation"),"label"] = 1-d.loc[d.split.eq("validation"),"label"]
        m=self.members[self.members.split.eq("train")]
        c=family_catalog(d,m)
        self.assertEqual(set(c.loc[c.eligible,"pfam"]),set(self.catalog.loc[self.catalog.eligible,"pfam"]))

    def test_ridge_svd_matches_sklearn(self):
        rng=np.random.default_rng(42)
        X=rng.normal(size=(38,9)); X[:,0]=3
        y=rng.normal(size=38)
        basis=RidgeBasis(X)
        for alpha in (.0001,1,10000):
            coef,intercept=basis.parameters(y,alpha)
            ref=Ridge(alpha=alpha,solver="svd").fit(basis.scaler.transform(X),y)
            np.testing.assert_allclose(coef,ref.coef_,atol=1e-12)
            np.testing.assert_allclose(intercept,ref.intercept_,atol=1e-12)

    def test_preprocessing_has_no_heldout_influence(self):
        rng=np.random.default_rng(123)
        X=rng.normal(size=(40,7)); y=rng.normal(size=40)
        fit=np.arange(40)<25
        c,i,s=fit_parameters(X,y,fit,1,False,{})
        X[~fit]=1e9; y[~fit]=-1e10
        c2,i2,s2=fit_parameters(X,y,fit,1,False,{})
        np.testing.assert_array_equal(s.mean_,s2.mean_)
        np.testing.assert_allclose(c,c2)
        self.assertEqual(i,i2)

    def test_weighted_regression_and_constant_target(self):
        y=np.array([1.,3.,7.,-2.]); p=np.array([2.,2.,4.,-1.])
        w=np.array([[1,2,3,1],[3,0,2,1]])
        result=regression_metrics(y,p,w)
        for j,weights in enumerate(w):
            self.assertAlmostEqual(result["r2"][j],r2_score(y,p,sample_weight=weights))
            self.assertAlmostEqual(result["mae"][j],mean_absolute_error(y,p,sample_weight=weights))
        self.assertTrue(np.isnan(regression_metrics(np.ones(4),p,w)["r2"]).all())
        self.assertTrue(np.isnan(binary_metrics(np.ones(4),p,w)["auroc"]).all())

    def test_small_fit_predictions_parameters_and_privacy(self):
        rng=np.random.default_rng(9)
        n=55
        X=rng.normal(size=(n,5)); y=X[:,0]*.6-X[:,2]*.2
        data=pd.DataFrame(dict(record_id=[f"fixture{i}" for i in range(n)],sha256=[str(i)*4 for i in range(n)],
                               split=["train"]*45+["validation"]*10,audit_cluster=np.arange(n),
                               inner_cv_fold=list(np.arange(45)%5)+[-1]*10))
        weights=np.ones((12,10),dtype=int)
        y[1]=np.nan
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)/"fit"
            fit_target(X,data,"fixture",False,y,folder,{},weights)
            self.assertTrue(valid_checkpoint(folder))
            table=pd.read_csv(folder/"validation_predictions.tsv",sep="\t")
            with np.load(folder/"parameters.npz") as params:
                pred=((X[45:]-params["mean"])/params["scale"])@params["coef"]+params["intercept"]
            np.testing.assert_allclose(table.prediction,pred)
            self.assertEqual(len(pd.read_csv(folder/"training_oof_predictions.tsv",sep="\t")),44)
            self.assertFalse(set(table.columns)&{"label","target_value","sequence","audit_cluster"})
            raw=make_bundle(folder)
            unpack_bundle(raw,Path(temp)/"restored")
            mark_local(Path(temp)/"restored",{"fixture":True})
            self.assertTrue(local_complete(Path(temp)/"restored",{"fixture":True}))
            (Path(temp)/"restored"/"parameters.npz").write_bytes(b"damaged")
            with self.assertRaises(AssertionError):
                local_complete(Path(temp)/"restored",{"fixture":True})

    def test_notebook_code_and_embedded_verification(self):
        import nbformat
        nb=build(); nbformat.validate(nb)
        for cell in nb.cells:
            if cell.cell_type != "code":
                continue
            compile(cell.source,"<cell>","exec")
            tree=ast.parse(cell.source)
            for node in ast.walk(tree):
                if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id=="run_stage1":
                    if node.args and isinstance(node.args[0],ast.Constant) and node.args[0].value=="-c":
                        compile(ast.literal_eval(node.args[1]),"<embedded-script>","exec")
        text="\n".join(c.source for c in nb.cells if c.cell_type=="code")
        self.assertNotIn("HR_TOKEN_PLM",text)
        self.assertNotIn("save_full_maps.py', '--",text)
        self.assertIn("data/pfam_memberships.tsv",text)

    def test_single_class_validation_is_retained_but_undefined(self):
        rng=np.random.default_rng(7)
        X=rng.normal(size=(60,3)); y=np.r_[np.arange(50)%2,np.ones(10)]
        d=pd.DataFrame(dict(record_id=[str(i) for i in range(60)],sha256=[str(i) for i in range(60)],
                            split=["train"]*50+["validation"]*10,audit_cluster=np.arange(60),
                            inner_cv_fold=list(np.arange(50)%5)+[-1]*10))
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp)
            fit_target(X,d,"single-class-fixture",True,y,folder,{},np.ones((10,10)))
            table=pd.read_csv(folder/"metrics.tsv",sep="\t")
            self.assertTrue(table.validation.isna().all())
            self.assertTrue(table.validation_n.eq(10).all())
            self.assertTrue(table.validation_positive.eq(10).all())
            self.assertEqual(len(pd.read_csv(folder/"validation_predictions.tsv",sep="\t")),10)

    def test_complete_report_smoke_and_population_guard(self):
        # Temporary, explicitly synthetic tables exercise every plot and aggregation path.
        from property_report import render
        from property_inputs import CONTINUOUS
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            cat=pd.DataFrame([dict(pfam="PF_fixture",name="Synthetic fixture",eligible=True,train_positive=8)])
            for layer in range(34):
                for target in ["solubility",*CONTINUOUS,"PF_fixture"]:
                    folder=root/f"layer_{layer:02d}"/target
                    folder.mkdir(parents=True)
                    keys=["r2","mae","rmse"] if target in CONTINUOUS else ["auroc","average_precision"]
                    rows=[dict(metric=k,validation=.6,low=.5,high=.7,cv_mean=.59,cv_sd=.03,
                               train_n=50,validation_n=10,validation_positive=5,parameter=1) for k in keys]
                    pd.DataFrame(rows).to_csv(folder/"metrics.tsv",sep="\t",index=False)
                    pd.DataFrame([dict(metric=k,validation=0,training_constant=.5) for k in keys]).to_csv(folder/"constant.tsv",sep="\t",index=False)
                    np.savez(folder/"bootstrap.npz",**{k:np.linspace(.5,.7,20) for k in keys})
            render(root,cat)
            self.assertEqual(len(list((root/"figures").glob("*.png"))),5)
            delta=pd.read_csv(root/"paired_vs_composition.tsv",sep="\t")
            self.assertTrue(delta.delta.eq(0).all())
            path=root/"layer_33"/"hydrophobicity"/"metrics.tsv"
            table=pd.read_csv(path,sep="\t"); table["validation_n"]=9
            table.to_csv(path,sep="\t",index=False)
            with self.assertRaises(AssertionError):
                render(root,cat)

    def test_completed_remote_cache_skips_archive_download(self):
        import hashlib
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); source=root/"source"; source.mkdir()
            (source/"result.tsv").write_text("record_id\tprediction\nfixture\t0.2\n")
            raw=make_bundle(source)
            state={"identity":{"fixture":1},"sha256":hashlib.sha256(raw).hexdigest()}
            store=PrivateStore.__new__(PrivateStore)
            store.names={"fixture/completed.json"}
            downloads=[]
            def get(name):
                downloads.append(name)
                return raw if name.endswith(".zip") else json.dumps(state).encode()
            store.get=get
            self.assertTrue(store.restore("fixture",root/"restored",{"fixture":1}))
            self.assertEqual(len(downloads),2)
            downloads.clear()
            self.assertTrue(store.restore("fixture",root/"restored",{"fixture":1}))
            self.assertEqual(downloads,["fixture/completed.json"])
            with self.assertRaises(AssertionError):
                store.restore("fixture",root/"restored",{"fixture":2})

    def test_existing_solubility_probe_is_reused_exactly(self):
        from probe import fit_probe
        rng=np.random.default_rng(5)
        X=rng.normal(size=(60,4))
        d=pd.DataFrame(dict(record_id=[str(i) for i in range(60)],sha256=[str(i) for i in range(60)],
                            split=["train"]*50+["validation"]*10,audit_cluster=np.arange(60),
                            inner_cv_fold=list(np.arange(50)%5)+[-1]*10,label=np.arange(60)%2))
        weights=np.ones((12,10),dtype=int)
        with tempfile.TemporaryDirectory() as temp:
            old=Path(temp)/"old"; old.mkdir()
            original=fit_probe(X,d,"composition",old,weights)
            new=Path(temp)/"new"
            reuse_solubility(old,"composition",d,new,weights)
            table=pd.read_csv(new/"metrics.tsv",sep="\t").set_index("metric")
            self.assertAlmostEqual(table.loc["auroc","validation"],original["auroc"])
            self.assertAlmostEqual(table.loc["auroc","cv_mean"],original["cv_mean"])
            self.assertEqual((old/"models/composition_parameters.npz").read_bytes(),(new/"parameters.npz").read_bytes())


if __name__ == "__main__":
    unittest.main()
