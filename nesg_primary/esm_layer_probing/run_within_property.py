"""Read completed global-probe predictions and evaluate subsets; never fit models."""
import argparse
import hashlib
import importlib.metadata as metadata
import io
import json
import os
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
from common import HERE,sha,write_json
from property_inputs import load_inputs
from run_property_probes import PrivateStore,mark_local
from run_layerwise_auroc import PINS,make_bundle
from within_property_evaluation import PROTOCOL,plan,read_predictions,evaluate_layer,render
from probe import bootstrap_weights

REQUIRED_INPUTS=["common.py","probe.py","property_inputs.py","property_probes.py","protocol.json",
                 "data/extraction_manifest.tsv","data/local_labels.tsv","data/input_integrity.json",
                 "data/pfam_memberships.tsv","data/uncontrolled/properties.tsv",
                 "data/uncontrolled/protocol.json","data/uncontrolled/integrity.json"]


class EvaluationStore(PrivateStore):
    def publish(self,prefix,folder,identity):
        from huggingface_hub import CommitOperationAdd
        from hf_storage import private_api
        private_api(self.repo)
        raw=make_bundle(folder)
        state=dict(identity=identity,sha256=hashlib.sha256(raw).hexdigest(),
                   observed_per_protein_targets_uploaded=False)
        self.remote(self.api.create_commit,repo_id=self.repo,repo_type="dataset",operations=[
            CommitOperationAdd(path_in_repo=prefix+"/results.zip",path_or_fileobj=raw),
            CommitOperationAdd(path_in_repo=prefix+"/completed.json",path_or_fileobj=json.dumps(state,sort_keys=True).encode()),
        ],commit_message="Save evaluation of unchanged NESG probes within biological strata")
        mark_local(folder,identity)
        print("Evaluation tables and figures saved to private HF.",flush=True)


def compatible_source(state):
    identity=state.get("identity",{})
    return (identity.get("train")==2361 and identity.get("validation")==691
            and all(identity.get("inputs",{}).get(name)==sha(HERE/name) for name in REQUIRED_INPUTS))


def select_source(store, requested=None):
    candidates=sorted(n.rsplit("/summary/completed.json",1)[0] for n in store.names
                      if n.startswith("esm_stage1/uncontrolled_property_probing_v1/") and n.endswith("/summary/completed.json"))
    if requested:
        requested=requested.rstrip("/")
        candidates=[p for p in candidates if p==requested]
    for prefix in candidates:
        state=json.loads(store.get(prefix+"/summary/completed.json"))
        if compatible_source(state):
            expected=[prefix+f"/layer_{i:02d}/completed.json" for i in range(34)]
            if all(n in store.names for n in expected):
                indices=[json.loads(store.get(n)) for n in expected]
                assert all(s["identity"]==state["identity"] for s in indices)
                return prefix,state["identity"],indices
    raise RuntimeError("No compatible complete uncontrolled-property run is available yet. Let that notebook finish through layer 33 and save its summary, then rerun this cell. No probes were fitted and no activations were downloaded. If you supplied SOURCE_RUN, check its exact HF prefix.")


def predictions_from_archive(raw, expected_sha, wanted, destination):
    """Verify the original checkpoint, extracting only selected prediction tables."""
    assert hashlib.sha256(raw).hexdigest()==expected_sha
    destination=Path(destination)
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        names=archive.namelist()
        assert len(names)==len(set(names))
        assert all(not Path(n).is_absolute() and ".." not in Path(n).parts for n in names)
        hashes=json.loads(archive.read("bundle_checksums.json"))
        assert set(names)==set(hashes)|{"bundle_checksums.json"}
        assert set(wanted)<=set(hashes)
        for name in wanted:
            content=archive.read(name)
            assert hashlib.sha256(content).hexdigest()==hashes[name]
            p=destination/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(content)


def prediction_folder(store,prefix,layer,index,source_identity,targets,cache,local_root):
    wanted=[name+"/validation_predictions.tsv" for name,_,_ in targets]
    digest=index["sha256"]
    local=local_root/prefix.rsplit("/",1)[1]/f"layer_{layer:02d}"
    marker=local/"local_complete.json"
    if marker.is_file():
        state=json.loads(marker.read_text())
        if state.get("identity")==source_identity and all(n in state.get("files",{}) and (local/n).is_file() and sha(local/n)==state["files"][n] for n in wanted):
            return local
    folder=cache/f"layer_{layer:02d}"
    marker=folder/"prediction_source.json"
    if marker.is_file():
        state=json.loads(marker.read_text())
        assert state["source_archive_sha256"]==digest
        assert set(state["files"])==set(wanted)
        assert all(sha(folder/n)==h for n,h in state["files"].items())
        return folder
    raw=store.get(prefix+f"/layer_{layer:02d}/results.zip")
    predictions_from_archive(raw,digest,wanted,folder)
    write_json(marker,dict(source_archive_sha256=digest,files={n:sha(folder/n) for n in wanted}))
    return folder


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",required=True)
    parser.add_argument("--source-run",default=None)
    parser.add_argument("--output",type=Path,default=Path("/content/nesg_within_property"))
    parser.add_argument("--existing-probe-output",type=Path,default=Path("/content/nesg_uncontrolled_probes"))
    args=parser.parse_args()
    data,properties,memberships,catalog=load_inputs()
    targets,jobs,support,edges,families=plan(data,properties,memberships,catalog)
    store=EvaluationStore(args.repo)
    source,source_identity,indices=select_source(store,args.source_run)
    identity=dict(protocol=PROTOCOL,source_run=source,source_identity=source_identity,
                  source_archive_sha256=[s["sha256"] for s in indices],
                  code_sha256={n:sha(HERE/n) for n in ["within_property_evaluation.py","run_within_property.py","run_property_probes.py","run_layerwise_auroc.py"]},
                  software={p:metadata.version(p) for p in ["numpy","pandas",*PINS]})
    fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    prefix=f"esm_stage1/within_property_evaluation_v1/{fingerprint}"
    result=args.output/fingerprint/"results"
    if store.restore(prefix,result,identity):
        write_json(args.output/"latest.json",dict(results=str(result),repository=args.repo,hf_prefix=prefix))
        print("Completed within-property evaluation reused; no metric recomputation.",flush=True)
        return
    result.mkdir(parents=True,exist_ok=True)
    support.to_csv(result/"stratum_target_support.tsv",sep="\t",index=False)
    families.to_csv(result/"all_training_family_support.tsv",sep="\t",index=False)
    write_json(result/"protocol.json",PROTOCOL)
    write_json(result/"training_quantile_edges.json",edges)
    write_json(result/"run_identity.json",identity)
    write_json(result/"source_provenance.json",dict(repository=args.repo,inspected_commit=store.info.sha,source_run=source))
    cache=args.output/fingerprint/"prediction_cache"
    def load(layer):
        folder=prediction_folder(store,source,layer,indices[layer],source_identity,targets,cache,args.existing_probe_output)
        return {name:read_predictions(folder/name/"validation_predictions.tsv",data,name,y) for name,_,y in targets}
    baseline=load(0)
    weights=bootstrap_weights(data.loc[data.split.eq("validation"),"audit_cluster"],PROTOCOL["bootstrap_repetitions"],PROTOCOL["seed"])
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=2):
        for layer in range(34):
            path=result/f"layer_{layer:02d}.tsv"; marker=result/f"layer_{layer:02d}_complete.json"
            if marker.is_file():
                state=json.loads(marker.read_text())
                assert state["fingerprint"]==fingerprint and sha(path)==state["sha256"]
                print(f"Reusing evaluated layer {layer:02d}",flush=True)
                continue
            scores=baseline if layer==0 else load(layer)
            table=evaluate_layer(layer,scores,baseline,jobs,weights)
            table.to_csv(path,sep="\t",index=False)
            write_json(marker,dict(fingerprint=fingerprint,sha256=sha(path)))
            print(f"Evaluated {'composition' if layer==0 else 'layer '+str(layer)}: saved predictions only; no fitting.",flush=True)
    render(result,support,catalog)
    write_json(result/"integrity.json",dict(layers=33,targets=len(targets),stratum_target_tasks=len(jobs),
               prediction_ids_hashes_splits_verified=True,refitted_models=0,ESM_forward_passes=0,reserved_test_accessed=False))
    store.publish(prefix,result,identity)
    write_json(args.output/"latest.json",dict(results=str(result),repository=args.repo,hf_prefix=prefix))
    print("Complete: all 33 layers evaluated within strata. See figures and support before interpreting changes.",flush=True)


if __name__=="__main__":
    os.environ.setdefault("MPLBACKEND","Agg")
    main()
