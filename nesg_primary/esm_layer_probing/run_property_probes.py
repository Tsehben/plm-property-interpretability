"""CPU-only uncontrolled probes, with verified per-layer private-HF checkpoints."""
import argparse
import hashlib
import importlib.metadata as metadata
import json
import os
from pathlib import Path
import shutil
import tempfile
from common import HERE, sha, token, write_json, composition
from run_layerwise_auroc import PINS, make_bundle, unpack_bundle
from property_inputs import load_inputs

CODE = ["common.py", "probe.py", "extract.py", "hf_storage.py", "save_full_maps.py",
        "run_layerwise_auroc.py", "property_inputs.py", "property_probes.py",
        "property_report.py", "run_property_probes.py"]


def inputs_identity():
    names = [*CODE, "protocol.json", "data/extraction_manifest.tsv", "data/local_labels.tsv",
             "data/pfam_memberships.tsv", "data/input_integrity.json"]
    names += [str(p.relative_to(HERE)) for p in sorted((HERE/"data/uncontrolled").iterdir()) if p.is_file()]
    return {name: sha(HERE/name) for name in names}


def local_complete(folder, identity):
    marker = folder/"local_complete.json"
    if not marker.is_file():
        return False
    state = json.loads(marker.read_text())
    assert state["identity"] == identity
    for name, digest in state["files"].items():
        assert not Path(name).is_absolute() and ".." not in Path(name).parts
        assert sha(folder/name) == digest, f"Damaged cached result: {name}"
    return True


def mark_local(folder, identity):
    write_json(folder/"local_complete.json",dict(identity=identity,files={
        str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob("*"))
        if p.is_file() and p.name != "local_complete.json"}))


def compatible_solubility(state, activation_hash):
    identity = state.get("identity", {})
    hashes = identity.get("inputs", {})
    required = ["common.py", "probe.py", "protocol.json", "data/extraction_manifest.tsv",
                "data/local_labels.tsv", "data/input_integrity.json"]
    return (state.get("models") == 34 and identity.get("activation_completion_sha256") == activation_hash
            and all(hashes.get(name) == sha(HERE/name) for name in required))


class PrivateStore:
    def __init__(self, repo):
        from hf_storage import private_api
        from save_full_maps import remote
        self.api, self.info = private_api(repo)
        self.repo, self.remote = repo, remote
        self.names = set(remote(self.api.list_repo_files, repo, repo_type="dataset", revision=self.info.sha))

    def get(self, name):
        from huggingface_hub import hf_hub_download
        return Path(self.remote(hf_hub_download,self.repo,name,repo_type="dataset",
                                revision=self.info.sha,token=token())).read_bytes()

    def restore(self, prefix, folder, identity):
        if prefix+"/completed.json" not in self.names:
            return False
        state = json.loads(self.get(prefix+"/completed.json"))
        assert state["identity"] == identity
        if local_complete(folder,identity):
            print(f"Reusing verified local {prefix.rsplit('/',1)[-1]}",flush=True)
            return True
        raw = self.get(prefix+"/results.zip")
        assert hashlib.sha256(raw).hexdigest() == state["sha256"]
        unpack_bundle(raw,folder)
        mark_local(folder,identity)
        print(f"Restored completed {prefix.rsplit('/',1)[-1]} from private HF",flush=True)
        return True

    def publish(self, prefix, folder, identity):
        from huggingface_hub import CommitOperationAdd
        raw = make_bundle(folder)
        state = dict(identity=identity,sha256=hashlib.sha256(raw).hexdigest(),
                     observed_per_protein_targets_uploaded=False)
        # No automatic repository creation. Privacy is rechecked before writing.
        from hf_storage import private_api
        private_api(self.repo)
        self.remote(self.api.create_commit,repo_id=self.repo,repo_type="dataset",operations=[
            CommitOperationAdd(path_in_repo=prefix+"/results.zip",path_or_fileobj=raw),
            CommitOperationAdd(path_in_repo=prefix+"/completed.json",path_or_fileobj=json.dumps(state,sort_keys=True).encode()),
        ],commit_message="Checkpoint NESG uncontrolled property probes")
        mark_local(folder,identity)
        print(f"Saved {prefix.rsplit('/',1)[-1]} to private HF",flush=True)

    def solubility(self, destination, activation_hash):
        candidates = sorted(n for n in self.names if n.startswith("esm_stage1/layerwise_probing_v1/") and n.endswith("/completed.json"))
        for name in candidates:
            state = json.loads(self.get(name))
            if not compatible_solubility(state,activation_hash):
                continue
            # Repository HEAD changes after our own checkpoints. Reuse depends on
            # the immutable bundle digest, not the current repository HEAD.
            old_identity = dict(repository=self.repo,index=name,sha256=state["sha256"])
            if not local_complete(destination,old_identity):
                raw = self.get(name.rsplit("/",1)[0]+"/results.zip")
                assert hashlib.sha256(raw).hexdigest() == state["sha256"]
                unpack_bundle(raw,destination)
                mark_local(destination,old_identity)
            print("Reusing 33 completed solubility probes and composition baseline; no refitting.",flush=True)
            return destination, old_identity
        print("No compatible completed solubility bundle found; solubility will use the same frozen protocol.",flush=True)
        return None, {"source": "new_fit_no_compatible_saved_run"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo")
    parser.add_argument("--activations",type=Path)
    parser.add_argument("--output",type=Path,default=Path("/content/nesg_uncontrolled_probes"))
    parser.add_argument("--threads",type=int,default=2)
    parser.add_argument("--baseline-only",action="store_true",help="Local actual composition fits; no HF/ESM access")
    args = parser.parse_args()
    if not args.baseline_only and not (args.repo and args.activations):
        parser.error("--repo and --activations required")
    assert args.threads > 0
    from threadpoolctl import threadpool_limits
    from property_probes import run_layer, target_list
    data, properties, memberships, catalog = load_inputs()
    targets = target_list(data,properties,memberships,catalog)
    identity = dict(inputs=inputs_identity(),software={p:metadata.version(p) for p in ["numpy","pandas",*PINS]},
                    threads=args.threads,train=2361,validation=691,reserved_test_accessed=False)
    if args.baseline_only:
        fingerprint = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
        folder = args.output/fingerprint/"layer_00"
        with threadpool_limits(limits=args.threads):
            run_layer(composition(data.sequence),data,targets,folder)
        write_json(args.output/"latest_baseline.json",dict(results=str(folder),identity=identity))
        return
    from hf_storage import verified_files, restore
    store = PrivateStore(args.repo)
    assert "esm_stage1/completed.json" in store.names, "Complete pooled extraction is required."
    # Read the small completion index first. A finished analysis needs no tensor download.
    identity["activation_completion_sha256"] = hashlib.sha256(store.get("esm_stage1/completed.json")).hexdigest()
    fingerprint = hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    prefix = f"esm_stage1/uncontrolled_property_probing_v1/{fingerprint}"
    root = args.output/fingerprint
    summary = root/"summary"
    if store.restore(prefix+"/summary",summary,identity):
        write_json(args.output/"latest.json",dict(results=str(summary),repository=args.repo,hf_prefix=prefix))
        print("All completed results reused. No models refitted.",flush=True)
        return
    old, reuse = store.solubility(root/"existing_solubility",identity["activation_completion_sha256"])
    write_json(root/"solubility_reuse.json",reuse)
    embeddings = None
    with threadpool_limits(limits=args.threads):
        for layer in range(34):
            folder = root/f"layer_{layer:02d}"
            if store.restore(prefix+f"/layer_{layer:02d}",folder,identity):
                continue
            if local_complete(folder,identity):
                store.publish(prefix+f"/layer_{layer:02d}",folder,identity)
                continue
            if layer and embeddings is None:
                from extract import load_all
                restore(args.repo,args.activations)
                verified_files(args.activations)
                assert sha(args.activations/"completed.json") == identity["activation_completion_sha256"]
                embeddings = load_all(args.activations,data)
                assert embeddings.shape == (3052,33,1280)
            X = composition(data.sequence) if layer == 0 else embeddings[:,layer-1,:]
            print(f"Fitting {'composition baseline' if layer == 0 else 'ESM layer '+str(layer)}: {len(targets)} targets",flush=True)
            run_layer(X,data,targets,folder,old,"composition" if layer == 0 else f"layer_{layer:02d}")
            store.publish(prefix+f"/layer_{layer:02d}",folder,identity)
    from property_report import render
    render(root,catalog)
    summary.mkdir(parents=True,exist_ok=True)
    for p in root.iterdir():
        if p.is_file() and p.suffix in {".tsv",".md",".json"}:
            shutil.copyfile(p,summary/p.name)
    shutil.copytree(root/"figures",summary/"figures",dirs_exist_ok=True)
    write_json(summary/"run_identity.json",identity)
    store.publish(prefix+"/summary",summary,identity)
    write_json(args.output/"latest.json",dict(results=str(summary),repository=args.repo,hf_prefix=prefix))
    print("Complete: uncontrolled probes only. Review the figures before controlled analyses.",flush=True)


if __name__ == "__main__":
    os.environ.setdefault("MPLBACKEND","Agg")
    main()
