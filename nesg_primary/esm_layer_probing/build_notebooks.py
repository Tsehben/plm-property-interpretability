"""Build portable Colab extraction and local, executable results notebooks."""

import inspect
import textwrap
import nbformat as nbf
from common import HERE


def md(text):
    return nbf.v4.new_markdown_cell(inspect.cleandoc(text))


def code(text):
    return nbf.v4.new_code_cell(textwrap.dedent(text).strip())


def notebook(cells):
    nb = nbf.v4.new_notebook(cells=cells)
    nb.metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python", "version": "3.12"},
    }
    return nb


def main(colab_only=False):
    colab = notebook(
        [
            md("""# NESG: frozen ESM-2 all-layer extraction
        **Reopening:** use Run all to check and resume. Installed packages, completed pooled activations and matching successful pilots are reused. Missing temporary files are restored from private HF; completed activations are not recomputed. Full-map extraction resumes committed batches when its optional cell is present.
        It extracts only 2,361 training + 691 validation proteins. Only the development extraction manifest is read; no reserved-test records or phenotype labels are loaded.
        Use a GPU runtime (T4 or better), enable notebook access to the Colab Secret **HF_TOKEN_PLM**, and mount the Google Drive containing your existing project folder in the first cell. No ZIP upload is needed. No model fine-tuning, probes or SAE analysis runs here.
        Hugging Face is the only persistent activation store. Every completed 64-protein shard is uploaded immediately to your explicitly configured private dataset. Google Drive supplies existing code and input sequences only. They are copied automatically to temporary Colab disk for execution; generated activations and model caches are never written to Drive.
        """),
            code("""
        from pathlib import Path
        import hashlib, shutil, os, subprocess, sys, json
        from google.colab import drive

        if not Path('/content/drive/MyDrive').is_dir():
            drive.mount('/content/drive')
        PROJECT_ROOT = Path('/content/drive/MyDrive/plm-property-interpretability')
        source = PROJECT_ROOT / 'nesg_primary/esm_layer_probing'
        work = Path('/content/nesg_esm_stage1')  # Temporary execution copy, not Drive.
        required = [
            'common.py', 'extract.py', 'hf_storage.py', 'requirements-colab.txt',
            'install_colab_dependencies.py', 'protocol.json',
            'data/extraction_manifest.tsv', 'data/input_integrity.json',
        ]
        missing = [name for name in required if not (source / name).is_file()]
        if missing:
            raise FileNotFoundError(f'Project files not found under {source}. Check PROJECT_ROOT and Drive sync. Missing: {missing}')
        for name in required:
            src = source / name
            dest = work / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            if not dest.is_file() or hashlib.sha256(dest.read_bytes()).digest() != hashlib.sha256(src.read_bytes()).digest():
                shutil.copyfile(src, dest)
            assert hashlib.sha256(dest.read_bytes()).digest() == hashlib.sha256(src.read_bytes()).digest()
        os.chdir(work)
        print('Read existing code and development inputs from:', source)
        print('Temporary execution directory:', work)
        print('Activation destination: private Hugging Face; no activation files will be saved to Drive.')
        """),
            md("""## 1. Check GPU, then prepare a separate extraction environment
        The installer first verifies a CUDA-enabled PyTorch build and an actual GPU allocation. A CPU-only build or disagreement between an actual import and its active metadata stops the run before any installation. Inactive duplicate metadata entries are logged and ignored; they do not force a runtime reset. It then creates `/content/nesg_esm_env`, which reads Colab's existing GPU/numerical libraries but keeps the pinned extraction packages separate from the notebook's packages. Only prebuilt wheels are allowed; no NumPy source builds. TensorFlow/Flax are disabled for this torch-only workflow.
        On reopening, pinned packages are checked in the dedicated interpreter. Matching installations skip pip; GPU/import checks still run. A deleted runtime needs a fresh environment, but saved HF results remain reusable."""),
            code("""
        import shlex
        from IPython.utils.process import system
        os.environ.update(USE_TORCH='1', USE_TF='0', USE_FLAX='0', PYTHONNOUSERSITE='1')
        status = system(shlex.join([sys.executable, '-u', 'install_colab_dependencies.py']))
        if status != 0:
            raise RuntimeError(f'Dependency setup stopped (exit {status}). See the live output above.')
        EXTRACTION_PYTHON = Path('/content/nesg_esm_env/bin/python')
        assert EXTRACTION_PYTHON.is_file()
        def run_stage1(*arguments):
            status = system(shlex.join([str(EXTRACTION_PYTHON), '-u', *map(str, arguments)]))
            if status != 0:
                raise RuntimeError(f'Stage 1 command stopped (exit {status}). See the live output above.')
        """),
            md("""## 2. Secure authentication and persistent output
        The token is retrieved only from Colab Secrets as `HF_TOKEN_PLM`, placed in the same environment variable, and passed directly to authenticated HF calls. It is never displayed or saved. Do not paste credentials into cells. The token needs write access to your configured private dataset in addition to model access. The suggested repository name is not proof that it exists. Leave `CREATE_PRIVATE_REPO=False` to check an existing private dataset; set it to `True` only to explicitly create the configured private dataset if missing. This requires organization creation permission. Public repositories are always refused."""),
            code("""
        from google.colab import userdata
        try:
            secret = userdata.get('HF_TOKEN_PLM')
        except Exception:
            raise RuntimeError('Enable notebook access to HF_TOKEN_PLM in Colab Secrets.') from None
        if not secret:
            raise RuntimeError('HF_TOKEN_PLM is missing from Colab Secrets.')
        os.environ['HF_TOKEN_PLM'] = secret
        del secret
        DESTINATION_REPO = 'PLM-Interp/nesg-esm2-layerwise'
        CREATE_PRIVATE_REPO = False  # Set True only to explicitly create this private dataset if missing.
        if not DESTINATION_REPO:
            raise ValueError('Set DESTINATION_REPO to your existing private HF dataset before continuing.')
        OUTPUT = Path('/content/nesg_esm_stage1_scratch/activations')
        OUTPUT.mkdir(parents=True, exist_ok=True)
        CACHE = '/content/nesg_esm_model_cache'  # Outside repository; can redownload after disconnection.
        # HF/model libraries run only in the dedicated interpreter, never the notebook kernel.
        if CREATE_PRIVATE_REPO:
            run_stage1('hf_storage.py', 'init', '--repo', DESTINATION_REPO, '--create-private')
        # Restore also checks authentication/privacy. It skips a verified complete
        # local cache and otherwise downloads only missing or damaged files.
        run_stage1('hf_storage.py', 'restore', '--repo', DESTINATION_REPO, '--folder', OUTPUT)
        """),
            md("""## 3. Real-model correctness pilot
        Fixed length quantiles plus all three X-containing proteins are chosen without labels.
        Checks cover token identities/counts, special-token exclusion, batching/padding invariance, 33 × 1,280 shapes, finite outputs and mixed precision against FP32 across every layer.
        Prespecified mixed tolerances: minimum cosine ≥0.999 and maximum relative L2 error ≤0.05. A failed pilot stops extraction; switch to FP32 rather than relaxing thresholds.
        FP32 padding check uses `atol=2e-4, rtol=2e-3`. A matching saved pilot is reused. If all pooled activations are already verified, no extraction pilot or model loading is needed."""),
            code("""
        PRECISION = 'float16'
        BATCH_SIZE = 2
        MAX_PADDED_TOKENS = 2048
        WORKERS = 0
        options = ['--output', str(OUTPUT), '--device', 'cuda', '--precision', PRECISION,
                   '--batch-size', str(BATCH_SIZE), '--max-tokens', str(MAX_PADDED_TOKENS),
                   '--workers', str(WORKERS), '--cache', CACHE, '--hf-repo', DESTINATION_REPO]
        run_stage1('extract.py', 'pilot', *options)
        if (OUTPUT / 'pilot.json').is_file():
            display(json.loads((OUTPUT / 'pilot.json').read_text()))
        else:
            print('Pooled extraction already complete; pilot rerun skipped.')
        """),
            md("""## 4. Extract all 33 transformer layers
        One forward pass per batch, then float32 residue means; full residue tensors are released immediately. `hidden_states[8]` means output at layer 8, not concatenated layers 1–8. Layer 0 is omitted (optional reference not requested for this run). Layer 33 follows the implementation's final encoder normalization. InterProt uses these same HF indices.
        If a session disconnects, rerun setup with the same private repository. Setup restores the uploaded shards into temporary storage and checks their hashes and protein alignment; extraction skips them. Each new shard is uploaded before the next batch of proteins starts. Upload failure stops the run; rerun after correcting connectivity or permissions. Only an unfinished or not-yet-uploaded shard may need recomputation. Reduce batch size on OOM; shard identities stay unchanged. Changing precision or protocol requires a different private repository and pilot; incompatible artifacts are never overwritten."""),
            code("""
        run_stage1('extract.py', 'extract', *options, '--shard-size', '64')
        """),
            md("""## 5. Validate complete output
        Expected 3,052 × 33 × 1,280 float32 values (~516 MB total, before metadata).
        The final completion manifest contains exact artifact checksums. Experimental labels and Pfam annotations stay in the separate local analysis manifest."""),
            code("""
        run_stage1('extract.py', 'validate', '--output', OUTPUT)
        completed = json.loads((OUTPUT / 'completed.json').read_text())
        print('Validated shape:', completed['shape'])
        print('Output:', OUTPUT)
        """),
            md("""## 6. Verify durable private storage
        Extraction uploads every validated shard immediately, then the final completion manifest. This cell rechecks privacy and retrieves the final manifest from the exact uploaded commit to confirm the run reached Hugging Face. Generated activations are not saved to Google Drive.
        Uploaded files contain activations and minimal ID/hash/split/extraction metadata, not sequences, experimental labels, Pfam annotations or reserved-test data."""),
            code(
                "verify_remote = \"import json, sys\\nfrom pathlib import Path\\nfrom common import token\\nfrom hf_storage import private_api\\nfrom huggingface_hub import hf_hub_download\\nrepo, output = sys.argv[1], Path(sys.argv[2])\\n_, repository = private_api(repo)\\ntry:\\n    remote_index = hf_hub_download(repo, 'esm_stage1/completed.json',\\n        repo_type='dataset', revision=repository.sha, token=token())\\nexcept Exception:\\n    raise RuntimeError('Could not verify remote completion; check token/repository/network. Details suppressed.') from None\\nassert json.loads(Path(remote_index).read_text()) == json.loads((output/'completed.json').read_text())\\nprint('All activations are stored in a verified private HF dataset:', repo)\\nprint('Remote commit:', repository.sha)\\n\"\nrun_stage1('-c', verify_remote, DESTINATION_REPO, OUTPUT)"
            ),
            md("""## 7. Reuse activations from Hugging Face
        The canonical activation copy is the private Hugging Face dataset. To run probes, use `hf_storage.py download --repo PLM-Interp/nesg-esm2-layerwise --folder /tmp/nesg_activations` with `HF_TOKEN_PLM`, then `probe.py --activations /tmp/nesg_activations`. Downloaded arrays are a disposable working cache, not an additional persistent backup.
        Existing experimental labels and Pfam annotations remain in the local research repository; they are not newly uploaded. Rerun the analysis notebook after probing. No model download is needed for probing. **Stop after Stage 1.**"""),
        ]
    )
    nbf.write(colab, HERE / "ESM_Extraction_Colab.ipynb")
    if colab_only:
        print(
            "Rebuilt Colab notebook: Drive inputs, HF activation storage; preserved executed analysis notebook."
        )
        return
    analysis = notebook(
        [
            md("""# NESG Stage 1: layer-wise linear solubility decoding
        This notebook separates completed measurements from pending ESM experiments. It reads only the frozen development cohort and saved results. No reserved-test sequences or labels, SAE features, interventions or biological-concept probes are accessed.
        The phenotype is high-expression NESG soluble recovery: E∈{4,5}, S∈{4,5} versus S=0. It is not a universal, condition-independent physical solubility measurement."""),
            code("""
        from pathlib import Path
        import sys, json, subprocess
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from IPython.display import display, Markdown
        here = Path.cwd()
        if not (here / 'common.py').exists():
            here = here / 'nesg_primary/esm_layer_probing'
        sys.path.insert(0, str(here))
        from common import HERE, PROTOCOL, read_local, sha
        from report import baseline_figure, layer_figure, recommendation, write_summary
        data = read_local()
        checks = json.loads((HERE / 'data/checks.json').read_text())
        display(pd.crosstab(data.split, data.label).rename(columns={0:'Insoluble',1:'Soluble'}))
        display(pd.Series(checks, name='Verified checks').to_frame())
        """),
            md("""## 1. Frozen inputs and protocol
        There are 2,361 training proteins (1,645 soluble / 716 insoluble) and 691 validation proteins (491 / 200). All three X-containing sequences remain unchanged; 50–810 residues fit without truncation. Stable IDs and SHA-256 values link the extraction manifest, labels, domain-level Pfam multi-memberships and MMseqs2 groups. Missing Pfam annotations never exclude proteins.
        Group assignment is homology-aware, not family-disjoint: a shared Pfam domain may appear in several splits. The prospective test lock remains; test labels were historically inspected in earlier reports and are not accessed here."""),
            code("""
        tr = data.split.eq('train'); va = data.split.eq('validation')
        assert set(data.loc[tr,'audit_cluster']).isdisjoint(data.loc[va,'audit_cluster'])
        assert data.loc[tr].groupby('audit_cluster').inner_cv_fold.nunique().max() == 1
        display(data.loc[tr].groupby('inner_cv_fold').agg(proteins=('record_id','size'),groups=('audit_cluster','nunique'),soluble=('label','sum')))
        display(PROTOCOL)
        """),
            md("""## 2. Fit or reuse the prespecified linear probes
        Each model is StandardScaler → L2 logistic regression, with no class balancing. C is selected by mean AUROC over the five frozen training folds, from [0.0001, 0.001, 0.01, 0.1, 1, 10]; exact ties prefer smaller C. Scaling is fitted inside every CV training fold. Selected-C fold scores are tuning estimates, not unbiased nested-CV estimates.
        All final models refit on the same 2,361 training proteins. Validation labels never fit or tune C. Set `REFIT=True` to regenerate saved results. Full ESM fitting runs only if a completed and validated extraction exists; absence of embeddings cannot generate a layer result."""),
            code("""
        REFIT = False
        activation_dir = HERE / 'activations'
        have_activations = (activation_dir / 'completed.json').exists()
        metrics_path = HERE / 'results/metrics.tsv'
        if REFIT or not metrics_path.exists():
            command = [sys.executable, str(HERE / 'probe.py')]
            if have_activations:
                command += ['--activations', str(activation_dir)]
            subprocess.run(command, check=True)
        metrics = pd.read_csv(metrics_path, sep='\\t')
        cv = pd.read_csv(HERE / 'results/all_cv.tsv', sep='\\t')
        status = json.loads((HERE / 'results/run_status.json').read_text())
        assert status['protocol_sha256'] == sha(HERE / 'protocol.json')
        assert status['input_manifest_sha256'] == sha(HERE / 'data/extraction_manifest.tsv')
        assert status['local_labels_sha256'] == sha(HERE / 'data/local_labels.tsv')
        full = (metrics.layer > 0).sum() == 33
        if have_activations and not full:
            display(Markdown('**New embeddings detected. Set REFIT=True or run probe.py --activations activations to fit the 33 ESM probes.**'))
        display(metrics)
        display(Markdown(f"**Completed: composition baseline; {int((metrics.layer > 0).sum())}/33 ESM probes.**"))
        """),
            md("""## 3. Sequence-composition diagnostic baseline
        Twenty standard amino-acid counts are divided by original sequence length, plus natural log length. X residues are retained in the denominator and not assigned to another amino acid; their unassigned mass is explicit. No aggregation or other biological-concept predictor enters the baseline.
        Confidence intervals use 2,000 resamples of the 364 validation homology groups, retaining all members of each sampled group. Bootstrap draws are shared by every model for paired comparisons and do not refit models."""),
            code("""
        baseline_figure(metrics, cv)
        plt.show()
        baseline = metrics[metrics.model.eq('composition')].iloc[0]
        display(Markdown(f"Composition validation AUROC: **{baseline.auroc:.4f}** (95% CI {baseline.auroc_low:.4f}–{baseline.auroc_high:.4f}); AP **{baseline.average_precision:.4f}**. Training-CV mean {baseline.cv_mean:.4f}, SD {baseline.cv_sd:.4f}. Validation soluble prevalence is {data.loc[va,'label'].mean():.4f}. Fold variation and bootstrap uncertainty limit interpretation of small performance changes."))
        """),
            md("""## 4. Figure 1 — Layer-wise validation AUROC
        Pointwise group-bootstrap intervals are conditional on the fitted probes and do not adjust for selecting the best of 33 layers. AUROC measures linear decodability under this protocol; it does not identify direct versus indirect encoding."""),
            code("""
        if full:
            layer_figure(1, metrics, cv, data.loc[va,'label'].mean()); plt.show()
            best = metrics[metrics.layer.gt(0)].sort_values('auroc', ascending=False).iloc[0]
            display(Markdown(f"Highest observed validation AUROC is layer {int(best.layer)}: {best.auroc:.4f} ({best.auroc_low:.4f}–{best.auroc_high:.4f}). This is a descriptive maximum across layers, not an independently confirmed optimum. Inspect paired differences and training-CV stability."))
        else:
            display(Markdown('**Not available:** the 33-layer ESM extraction/probing has not run. No layer curve is plotted.'))
        """),
            md("""## 5. Figure 2 — Training-CV consistency
        Every layer uses identical proteins and group folds. Shared fold trajectories can reflect biological differences among folds as well as sampling variability. The displayed scores reuse hyperparameter-selection folds."""),
            code("""
        if full:
            layer_figure(2, metrics, cv, data.loc[va,'label'].mean()); plt.show()
            display(metrics[metrics.layer.gt(0)][['layer','C','cv_mean','cv_sd','cv_min','cv_max']])
            display(Markdown('Compare depth trends across folds before treating a validation peak as reproducible. CV mean and dispersion are both retained for later layer selection.'))
        else:
            display(Markdown('**Pending ESM results.** Composition CV estimates are available above.'))
        """),
            md("""## 6. Figure 3 — Layer-wise average precision
        Average precision depends on class prevalence. The reference is the same 491/691 validation soluble prevalence for all models; validation is never balanced."""),
            code("""
        if full:
            layer_figure(3, metrics, cv, data.loc[va,'label'].mean()); plt.show()
            display(Markdown('AP complements AUROC for this naturally imbalanced cohort. It is not comparable across populations with different prevalence without further qualification.'))
        else:
            display(Markdown('**Pending ESM results.** No AP curve is plotted.'))
        """),
            md("""## 7. Figure 4 — Prespecified representative layers
        Layers 8, 16 and 32 span early, middle and late depth and have released InterProt checkpoints. They were fixed before ESM results. Paired intervals use identical group-bootstrap resamples for both layers. An interval crossing zero indicates the observed difference is not clearly resolved under this resampling scheme, not proof of equivalence."""),
            code("""
        if full:
            layer_figure(4, metrics, cv, data.loc[va,'label'].mean()); plt.show()
            display(pd.read_csv(HERE / 'results/representative_paired_differences.tsv', sep='\\t'))
            display(metrics[metrics.layer.gt(0)][['layer','delta_vs_composition','delta_low','delta_high']])
        else:
            display(Markdown('**Pending ESM results.** Representative-layer comparisons and paired baseline differences cannot yet be estimated.'))
        """),
            md("""## 8. Integrity, interpretation and Stage 2 handoff
        Verify that all actual predictions use the same validation ID/hash/label/group rows. Multi-Pfam memberships and overlapping-domain flags remain in separate local tables; no within-family AUROC is calculated here.
        Candidate selection after completion considers one available InterProt layer per depth bin (1–11, 12–22, 23–33). Within each bin, retain candidates within one fold SD of the highest training-CV mean, then prefer lower fold SD (ties: earlier layer). This is a pragmatic depth/stability rule, not a statistical confidence set; it does not select the top three validation AUROCs."""),
            code("""
        expected = data.loc[va, ['record_id','sha256','split','label','audit_cluster','inner_cv_fold']].reset_index(drop=True)
        for model in metrics.model:
            pred = pd.read_csv(HERE / f'results/{model}_validation_predictions.tsv', sep='\\t')
            pd.testing.assert_frame_equal(pred[expected.columns], expected, check_dtype=False)
            assert np.isfinite(pred[['decision_score','probability']]).all().all()
        if full:
            display(Markdown(f'Candidate SAE layers for a subsequent, separately authorized stage: **{recommendation(metrics)}**. These are depth-spanning candidates, not proof of mechanism.'))
        else:
            display(Markdown('**No evidence-based SAE-layer recommendation yet.** 8/16/32 are only prespecified comparison anchors. Complete Colab extraction and the 33 probes first.'))
        write_summary(metrics, cv)
        print('Prediction alignment checks passed. No reserved-test evaluation. Stage 2 not run.')
        """),
            md("""## Reproduction and references
        Model and tokenizer: [pinned ESM-2 revision](https://huggingface.co/facebook/esm2_t33_650M_UR50D/commit/08e4846e537177426273712802403f7ba8261b6c).
        [InterProt model card](https://huggingface.co/liambai/InterProt-ESM2-SAEs/blob/3ded2e51641df84d6b89f25430874e3a8eb42c24/README.md) uses the same `hidden_states[LAYER]` convention. Layer 33 includes final encoder normalization; intermediate outputs are not concatenations.
        Dependencies are pinned in `requirements.txt`; executed environment and artifact hashes are recorded by `execute_analysis.py`. See `README.md` for commands and private storage safeguards. Current partial execution does not imply the real-model extraction pilot has passed."""),
        ]
    )
    nbf.write(analysis, HERE / "ESM_Layerwise_Solubility.ipynb")
    print("Built Colab and analysis notebooks.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--colab-only", action="store_true")
    main(colab_only=parser.parse_args().colab_only)
