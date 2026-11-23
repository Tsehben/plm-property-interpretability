"""Evaluate saved global probes within strata. No model fitting or ESM inference."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from common import write_json
from property_inputs import CONTINUOUS
from property_probes import target_list, scalar_metrics, binary_metrics, regression_metrics
from probe import bootstrap_weights

PROTOCOL = dict(version=1, layers=list(range(1,34)), evaluation_split="validation",
                scope="unchanged global probes evaluated within one restriction at a time",
                numeric_quantiles=[.2,.4,.6,.8], quantile_method="linear; observed training proteins only",
                interval_rule="lower inclusive, upper exclusive; outer intervals unbounded; tied edges collapsed",
                family_definitions="all GA-passing memberships; separate noncompeting-assignment sensitivity",
                family_catalog="families observed in training; no training solubility-support requirement",
                family_point_metrics_min_validation_n=2, bootstrap_repetitions=2000, seed=20261118,
                bootstrap_scope="solubility only; same validation-group draws across layers and baseline",
                min_groups_per_solubility_class_for_ci=3, ci_min_valid_fraction=.9,
                secondary_targets="continuous targets within families; all existing targets within numeric ranges",
                reserved_test_accessed=False, refit=False, layer_selection=False)


def make_strata(data, properties, memberships):
    tr=data.split.eq("train").to_numpy(); va=~tr
    strata=[dict(key="overall",kind="overall",control="none",category="All validation",variant="all",mask=np.ones(len(data),bool))]
    edges={}
    for control in CONTINUOUS:
        values=properties[control].to_numpy(float)
        observed=values[tr & np.isfinite(values)]
        cuts=np.unique(np.quantile(observed,PROTOCOL["numeric_quantiles"],method="linear"))
        cuts=cuts[(cuts>observed.min()) & (cuts<observed.max())]
        edges[control]=cuts.tolist()
        assignment=np.searchsorted(cuts,values,side="right")
        for i in range(len(cuts)+1):
            mask=np.isfinite(values) & (assignment==i)
            strata.append(dict(key=f"range:{control}:{i+1}",kind="range",control=control,
                               category=f"Range {i+1}",variant="training_quantiles",mask=mask,
                               lower=float(cuts[i-1]) if i else -np.inf,
                               upper=float(cuts[i]) if i<len(cuts) else np.inf))
    family_rows=[]
    train_ids=set(data.loc[tr,"record_id"])
    for family, group in memberships.groupby("pfam",sort=True):
        if not set(group.record_id)&train_ids:
            continue
        first=group.iloc[0]
        flagged=group.competing_overlap.astype(str).str.lower().eq("true")
        for variant, g in [("all_ga",group),("noncompeting",group[~flagged])]:
            mask=data.record_id.isin(g.record_id).to_numpy()
            record=dict(key=f"family:{family}:{variant}",kind="family",control=family,
                        category=str(first["name"]),variant=variant,mask=mask)
            family_rows.append(dict(pfam=family,name=first["name"],variant=variant,
                                    train_n=int((tr&mask).sum()),validation_n=int((va&mask).sum()),
                                    validation_soluble=int((va&mask&data.label.eq(1).to_numpy()).sum()),
                                    validation_insoluble=int((va&mask&data.label.eq(0).to_numpy()).sum())))
            # Singleton/absent families remain in the support table; no discrimination can be estimated.
            if (va&mask).sum() >= PROTOCOL["family_point_metrics_min_validation_n"]:
                strata.append(record)
    return strata,edges,pd.DataFrame(family_rows)


def plan(data, properties, memberships, catalog):
    strata,edges,families=make_strata(data,properties,memberships)
    targets=target_list(data,properties,memberships,catalog)
    tr=data.split.eq("train").to_numpy(); va=~tr
    jobs=[]; support=[]
    for target,binary,y in targets:
        finite=np.isfinite(y)
        for stratum in strata:
            # No family-membership task inside a defining family; family targets are evaluated in numeric ranges only.
            if stratum["kind"]=="family" and target not in ["solubility",*CONTINUOUS]:
                continue
            tmask=tr & finite & stratum["mask"]
            vmask=va & finite & stratum["mask"]
            n=int(vmask.sum()); outcome=y[vmask]
            groups=data.loc[vmask,"audit_cluster"].to_numpy()
            positive_groups=len(np.unique(groups[outcome==1])) if binary else None
            negative_groups=len(np.unique(groups[outcome==0])) if binary else None
            reason="estimable"
            if n<2:
                reason="fewer_than_two_observed_validation_proteins"
            elif binary and np.unique(outcome).size<2:
                reason="single_validation_class"
            elif not binary and np.ptp(outcome)==0:
                reason="constant_validation_target_R2_undefined"
            row={k:v for k,v in stratum.items() if k!="mask"}
            row.update(target=target,binary=binary,train_n=int(tmask.sum()),validation_n=n,
                       validation_groups=len(np.unique(groups)),validation_positive=int(outcome.sum()) if binary else np.nan,
                       positive_groups=positive_groups,negative_groups=negative_groups,
                       validation_target_sd=float(np.std(outcome)) if n else np.nan,
                       validation_target_min=float(outcome.min()) if n else np.nan,
                       validation_target_max=float(outcome.max()) if n else np.nan,status=reason,
                       same_property_range=(stratum["control"]==target))
            ci=(target=="solubility" and reason=="estimable" and min(positive_groups,negative_groups)>=PROTOCOL["min_groups_per_solubility_class_for_ci"])
            row["ci_support_eligible"]=ci
            jobs.append(dict(info=row,mask=vmask,within_validation=vmask[va],outcome=outcome,binary=binary))
            support.append(row)
    return targets,jobs,pd.DataFrame(support),edges,families


def read_predictions(path,data,target,y):
    """Reject missing/extra/duplicate or hash-mismatched records before any metric calculation."""
    expected=data.loc[data.split.eq("validation") & np.isfinite(y),["record_id","sha256","split"]]
    table=pd.read_csv(path,sep="\t",dtype={"record_id":str,"sha256":str,"split":str},float_precision="round_trip")
    assert table.record_id.is_unique and set(table.record_id)==set(expected.record_id), f"Prediction IDs differ: {target}"
    table=table.set_index("record_id").loc[expected.record_id]
    assert table.sha256.tolist()==expected.sha256.tolist(), f"Prediction hashes differ: {target}"
    assert table.split.eq("validation").all()
    assert np.isfinite(table.prediction.to_numpy()).all()
    out=np.full(len(data),np.nan)
    out[expected.index]=table.prediction.to_numpy()
    return out


def evaluate_layer(layer, predictions, baseline, jobs, weights):
    rows=[]
    for job in jobs:
        info=job["info"]; target=info["target"]; mask=job["mask"]
        y=job["outcome"]; p=predictions[target][mask]; b=baseline[target][mask]
        assert np.isfinite(p).all() and np.isfinite(b).all()
        keys=["auroc","average_precision"] if job["binary"] else ["r2","mae","rmse"]
        values=scalar_metrics(y,p,job["binary"]) if len(y)>=2 else {k:np.nan for k in keys}
        reference=scalar_metrics(y,b,job["binary"]) if len(y)>=2 else {k:np.nan for k in keys}
        boot={}; bboot={}
        if info["ci_support_eligible"]:
            w=weights[:,job["within_validation"]]
            boot=binary_metrics(y,p,w); bboot=binary_metrics(y,b,w)
        for metric in keys:
            low=high=dlo=dhi=np.nan; valid=0
            if metric in boot:
                good=np.isfinite(boot[metric]) & np.isfinite(bboot[metric])
                valid=int(good.sum())
                if good.mean()>=PROTOCOL["ci_min_valid_fraction"]:
                    low,high=np.quantile(boot[metric][good],[.025,.975])
                    dlo,dhi=np.quantile((boot[metric]-bboot[metric])[good],[.025,.975])
            rows.append(dict(key=info["key"],target=target,layer=layer,metric=metric,value=values[metric],
                             low=low,high=high,composition=reference[metric],
                             delta_vs_composition=values[metric]-reference[metric],delta_low=dlo,delta_high=dhi,
                             bootstrap_valid=valid))
    return pd.DataFrame(rows)


def render(destination, support, catalog):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    destination=Path(destination)
    metrics=pd.concat([pd.read_csv(destination/f"layer_{i:02d}.tsv",sep="\t") for i in range(34)],ignore_index=True)
    assert set(metrics.layer)==set(range(34))
    table=metrics.merge(support,on=["key","target"],validate="many_to_one")
    assert len(table)==len(metrics)
    table.to_csv(destination/"within_property_metrics.tsv",sep="\t",index=False)
    figures=destination/"figures"; figures.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    def save(fig,name):
        fig.tight_layout()
        for ext in ("png","svg"):
            fig.savefig(figures/f"{name}.{ext}",dpi=160,bbox_inches="tight")
        plt.close(fig)
    sol=table[(table.target=="solubility") & (table.metric=="auroc") & table.layer.gt(0)]
    fig,axes=plt.subplots(3,2,figsize=(13,12))
    full=sol[sol.key.eq("overall")].sort_values("layer")
    for ax,control in zip(axes.flat,CONTINUOUS):
        ax.plot(full.layer,full.value,"k--",label="All validation (different population)",alpha=.65)
        for key,g in sol[(sol.kind=="range") & sol.control.eq(control)].groupby("key",sort=True):
            g=g.sort_values("layer"); row=g.iloc[0]
            line=ax.plot(g.layer,g.value,label=f"{row.category}: n={int(row.validation_n)}; +{int(row.validation_positive)}/−{int(row.validation_n-row.validation_positive)}")[0]
            ax.fill_between(g.layer,g.low,g.high,color=line.get_color(),alpha=.08)
        ax.axhline(.5,color="gray",ls=":")
        ax.set(title=control,xlabel="ESM-2 layer",ylabel="Validation solubility AUROC",ylim=(0,1),xticks=[1,8,16,24,33])
        ax.legend(fontsize=7)
    axes.flat[-1].axis("off")
    axes.flat[-1].text(0,.9,"Same globally trained probes; no refitting.\nOne property restricted at a time.\nRanges use training-only quintiles.\nBands: conditional group-bootstrap 95% CI.\nRange restriction is not exact matching.\nAcross-population differences are descriptive.",va="top",linespacing=1.8)
    save(fig,"solubility_within_ranges")
    for variant in ("all_ga","noncompeting"):
        sub=sol[sol.kind.eq("family") & sol.variant.eq(variant)]
        sub=sub[np.isfinite(sub.value)]
        fig,ax=plt.subplots(figsize=(13,max(4,sub.key.nunique()*.29)))
        if len(sub):
            descriptions=sub.drop_duplicates("key").sort_values(["train_n","control"],ascending=[False,True]).set_index("key")
            mat=sub.pivot(index="key",columns="layer",values="value").reindex(descriptions.index)
            labels=[f"{r.control} {r.category}: +{int(r.validation_positive)}/−{int(r.validation_n-r.validation_positive)}; groups {int(r.positive_groups)}/{int(r.negative_groups)}" for r in descriptions.itertuples()]
            im=ax.imshow(np.ma.masked_invalid(mat.to_numpy()),aspect="auto",vmin=0,vmax=1,cmap="viridis")
            ax.set(yticks=range(len(mat)),yticklabels=labels,xticks=[0,7,15,23,32],xticklabels=[1,8,16,24,33],xlabel="ESM-2 layer")
            fig.colorbar(im,ax=ax,label="Validation AUROC")
        else:
            ax.text(.5,.5,"No estimable mixed-label families",ha="center",transform=ax.transAxes)
        ax.set_title(f"Within-family soluble recovery — {variant}\nCounts and groups are soluble/insoluble. Sparse rows are exploratory; overlapping families are not independent.")
        save(fig,f"solubility_within_families_{variant}")
    # Secondary property heatmaps include every training-defined numeric stratum.
    for target in CONTINUOUS:
        sub=table[(table.target==target) & (table.metric=="r2") & table.kind.eq("range") & table.layer.gt(0)]
        pivot=sub.pivot(index="key",columns="layer",values="value").sort_index()
        fig,ax=plt.subplots(figsize=(13,9))
        im=ax.imshow(np.ma.masked_invalid(pivot.to_numpy()),aspect="auto",vmin=-1,vmax=1,cmap="coolwarm_r")
        ax.set(yticks=range(len(pivot)),yticklabels=[k.replace("range:","") for k in pivot.index],
               xticks=[0,7,15,23,32],xticklabels=[1,8,16,24,33],xlabel="ESM-2 layer",
               title=f"{target}: validation R² within property ranges\nColor clips below −1; exact scores and target SD remain in tables. No regression CI.")
        fig.colorbar(im,ax=ax,label="R² (higher is better)")
        save(fig,f"secondary_{target}")
    (destination/"WITHIN_PROPERTY_RESULTS.md").write_text("""# Within-property evaluation

These results evaluate the same globally trained probes at all 33 layers on fixed validation subsets. There was no fitting, ESM inference, layer selection or reserved-test access. Layer 00 is the existing composition baseline, not the ESM embedding layer.

Numeric strata use observed training quintiles with tied boundaries collapsed and unbounded outer intervals. They restrict one descriptor at a time; they do not hold that descriptor exactly constant or control all other properties. Missing descriptors remain missing. Family strata retain multiple GA-passing memberships, with noncompeting-assignment results reported separately on their own populations. Domain membership does not guarantee matched architecture, organism or assay.

Every family observed in training is included in the support table. Families with at least two validation proteins receive point estimates where the target varies. Single-class binary targets have undefined AUROC/AP; constant continuous targets have undefined R². Small families may give extreme AUROC values by chance. These strata are not new training datasets. The defining family is not predicted within itself.

Solubility intervals resample validation homology groups and are conditional on the saved probes. They require ≥3 groups per label and ≥90% valid draws; these are reporting safeguards, not proof of adequate confirmatory power. Secondary targets have point estimates only. Exact target range and SD accompany regression results because range restriction changes R². Existing family-membership probes are assessed within numeric ranges only; other properties are also assessed within families.

`within_property_metrics.tsv` contains all targets, layers, strata, support and metrics. `delta_vs_composition` and its interval compare models on identical proteins within a stratum. Overall-versus-stratum curves use different populations and are descriptive, not paired treatment effects or causal evidence of confounding. No families are pooled into a duplicated protein population. Plotting many strata/layers is exploratory; the largest drop is not a prespecified finding.

Training counts describe support. Performance here is fixed-validation performance; training probes are unchanged. Review uncontrolled results first, then these evaluations. Choosing a few layers may simplify a slide, but all 33 remain the primary analysis.
""")
