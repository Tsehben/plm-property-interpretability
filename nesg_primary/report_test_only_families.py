"""User-authorized test-family inspection using label-free split IDs and cached Pfam.

No solubility labels, sequences, model predictions, or restricted test manifests
are opened. This inspects held-out biological composition, not model performance.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import pandas as pd

ROOT=Path(__file__).resolve().parent


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def summarize(splits,domains):
    members=domains[["sha256","pfam"]].drop_duplicates().merge(
        splits[["record_id","sha256","audit_cluster","split"]],on="sha256",validate="many_to_one")
    assert not members.duplicated(["record_id","pfam"]).any()
    sets={s:set(members.loc[members.split.eq(s),"pfam"]) for s in ["train","validation","reserved_test"]}
    only=sets["reserved_test"]-(sets["train"]|sets["validation"])
    test=members[members.split.eq("reserved_test")]
    novel=test[test.pfam.isin(only)]
    table=novel.groupby("pfam").agg(test_unique_proteins=("record_id","nunique"),
                                    test_homology_groups=("audit_cluster","nunique")).reset_index()
    meta=domains[["pfam","family_name","description","clan"]].drop_duplicates()
    assert meta.pfam.is_unique
    table=table.merge(meta,on="pfam",validate="one_to_one")
    evidence=domains.groupby(["sha256","pfam"]).competing_overlap.any().reset_index()
    flagged=novel.merge(evidence,on=["sha256","pfam"],validate="one_to_one")
    counts=flagged[flagged.competing_overlap].groupby("pfam").record_id.nunique()
    table["test_proteins_with_competing_assignments"]=table.pfam.map(counts).fillna(0).astype(int)
    table["train_unique_proteins"]=0; table["validation_unique_proteins"]=0
    table=table.sort_values(["test_unique_proteins","pfam"],ascending=[False,True]).reset_index(drop=True)
    table=table[["pfam","family_name","description","clan","test_unique_proteins","test_homology_groups",
                 "test_proteins_with_competing_assignments","train_unique_proteins","validation_unique_proteins"]]
    test_n=int(splits.split.eq("reserved_test").sum())
    novel_ids=set(novel.record_id)
    known_ids=set(test.loc[~test.pfam.isin(only),"record_id"])
    annotated=set(test.record_id)
    stats=dict(distinct_families_by_split={s:len(v) for s,v in sets.items()},
               total_cohort_distinct_families=len(set.union(*sets.values())),
               test_families_absent_from_train=len(sets["reserved_test"]-sets["train"]),
               test_families_absent_from_train_but_present_in_validation=len((sets["reserved_test"]-sets["train"])&sets["validation"]),
               test_only_families=len(only),reserved_test_proteins=test_n,
               test_only_families_pct_of_test_families=100*len(only)/len(sets["reserved_test"]),
               test_proteins_with_any_test_only_family=len(novel_ids),
               test_proteins_with_any_test_only_family_pct=100*len(novel_ids)/test_n,
               test_proteins_with_exclusively_test_only_families=len(novel_ids-known_ids),
               test_proteins_with_both_test_only_and_development_families=len(novel_ids&known_ids),
               test_proteins_with_only_development_families=len(annotated-novel_ids),
               test_proteins_without_pfam=test_n-len(annotated))
    assert len(table)==len(only)
    assert (len(novel_ids-known_ids)+len(novel_ids&known_ids)+len(annotated-novel_ids)+test_n-len(annotated))==test_n
    return table,stats


def main():
    split_path=ROOT/"results/split_assignments.tsv"
    domains_path=ROOT.parent/"dataset_comparison/results/domains_with_clans.tsv.gz"
    lock=json.loads((ROOT/"results/test_lock.json").read_text())
    pins=json.loads((ROOT.parent/"dataset_comparison/results/output_checksums.json").read_text())
    before={str(p.relative_to(ROOT.parent)):digest(p) for p in [split_path,domains_path]}
    assert digest(split_path)==lock["rerun_safe_checksums"]["split_assignments.tsv"]
    assert digest(domains_path)==pins["results/domains_with_clans.tsv.gz"]
    splits=pd.read_csv(split_path,sep="\t",usecols=["record_id","sha256","audit_cluster","split"])
    assert splits.record_id.is_unique and splits.sha256.is_unique
    assert splits.split.value_counts().to_dict()=={"train":2361,"reserved_test":710,"validation":691}
    cols=["sha256","pfam","family_name","description","clan","competing_overlap",
          "sequence_score","domain_score","sequence_ga","domain_ga"]
    domains=pd.read_csv(domains_path,sep="\t",usecols=cols,keep_default_na=False)
    domains=domains[domains.sha256.isin(splits.sha256)].copy()
    assert domains.competing_overlap.isin([True,False]).all()
    assert domains.sequence_score.ge(domains.sequence_ga).all() and domains.domain_score.ge(domains.domain_ga).all()
    table,stats=summarize(splits,domains)
    assert stats["distinct_families_by_split"]=={"train":949,"validation":402,"reserved_test":418}
    assert stats["total_cohort_distinct_families"]==1383
    out=ROOT/"results/test_only_families"; out.mkdir(parents=True,exist_ok=True)
    table.to_csv(out/"test_only_families.tsv",sep="\t",index=False)
    table.to_csv(out/"test_only_families.csv",index=False)
    (out/"summary.json").write_text(json.dumps(stats,indent=2)+"\n")
    assert before=={str(p.relative_to(ROOT.parent)):digest(p) for p in [split_path,domains_path]}
    event=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),
               authorization="User requested checking whether any families occur only in the test set",
               scope="held-out family composition inspection, all GA-passing Pfam memberships",
               solubility_labels_loaded=False,sequences_loaded=False,model_predictions_loaded=False,
               restricted_test_manifest_opened=False,UESolDS_independent_test_opened=False,
               cohort_splits_annotations_unchanged=True,source_sha256=before,
               output_sha256={p.name:digest(p) for p in out.iterdir() if p.name!="inspection_event.json"},
               interpretation="Descriptive inspection only; do not use to revise splits, model selection or hypotheses")
    (out/"inspection_event.json").write_text(json.dumps(event,indent=2)+"\n")
    print(json.dumps(stats,indent=2))
    print(table.head(12).to_string(index=False))


if __name__=="__main__": main()
