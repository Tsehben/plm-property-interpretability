"""Build the standalone Colab notebook without changing extraction/analysis notebooks."""

from pathlib import Path
import nbformat
from build_notebooks import code, md, notebook


def build():
    nb = notebook(
        [
            md("""# NESG: train and validation AUROC across ESM-2 layers

        Open this notebook in Colab and choose **Run all**. A CPU runtime is sufficient;
        no ESM or SAE inference is performed. Existing pooled activations are restored
        from the private HF dataset. Full residue maps are not downloaded.

        The notebook fits 33 standardized L2 logistic probes and a composition baseline,
        using the frozen 2,361 training / 691 validation proteins and training-only
        homology-aware folds. It displays both train-CV versus validation and train-fit
        versus validation figures. The reserved test is not accessed.

        **Reopening:** installed packages and valid cached files are reused. Complete
        matching results on HF are downloaded instead of refitted. Partial fits resume
        within the same runtime; deleting a runtime before final upload loses those
        partial fits, but the saved ESM activations remain available.

        **Execution status:** the workflow has local engineering tests; the actual 33-layer
        results will appear when run with your Colab Secret `HF_TOKEN_PLM`. No placeholder
        AUROC curves are included. HF remains the persistent result store; Drive supplies
        code and the frozen local development manifest."""),
            md("""## 1. Locate the project and copy only changed inputs

        Change `PROJECT_ROOT` only if your Drive project is in a different location.
        No credential files, full-map tensors or reserved-test files are copied."""),
            code("""
        from pathlib import Path
        import hashlib, json, os, shutil, subprocess, sys, venv
        from google.colab import drive

        if not Path('/content/drive/MyDrive').is_dir():
            drive.mount('/content/drive')
        PROJECT_ROOT = Path('/content/drive/MyDrive/plm-property-interpretability')
        source = PROJECT_ROOT / 'nesg_primary/esm_layer_probing'
        work = Path('/content/nesg_esm_stage1')
        required = [
            'common.py', 'extract.py', 'hf_storage.py', 'save_full_maps.py',
            'probe.py', 'plot_train_validation.py', 'run_layerwise_auroc.py',
            'protocol.json', 'data/extraction_manifest.tsv',
            'data/input_integrity.json', 'data/local_labels.tsv',
        ]
        missing = [name for name in required if not (source / name).is_file()]
        if missing:
            raise FileNotFoundError(f'Check PROJECT_ROOT and Drive sync. Missing: {missing}')
        changed = 0
        for name in required:
            src, dst = source / name, work / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            digest = hashlib.sha256(src.read_bytes()).digest()
            if not dst.is_file() or hashlib.sha256(dst.read_bytes()).digest() != digest:
                shutil.copyfile(src, dst)
                changed += 1
            assert hashlib.sha256(dst.read_bytes()).digest() == digest
        os.chdir(work)
        print(f'Updated {changed} code/input files; unchanged files reused.')
        """),
            md("""## 2. Prepare or reuse the CPU probing environment

        Reuse the extraction interpreter if it still exists. Otherwise create a small
        separate environment inheriting Colab's numerical stack. Only missing probing
        packages are installed as prebuilt wheels; there is no CUDA preflight, ESM
        download or PyTorch installation. Package versions are recorded in the results."""),
            code("""
        import shlex
        from IPython.utils.process import system

        os.environ.update(USE_TORCH='1', USE_TF='0', USE_FLAX='0', PYTHONNOUSERSITE='1', MPLBACKEND='Agg')
        EXTRACTION_PYTHON = Path('/content/nesg_esm_env/bin/python')
        if not EXTRACTION_PYTHON.is_file():
            environment = Path('/content/nesg_probe_env')
            EXTRACTION_PYTHON = environment / 'bin/python'
            if not EXTRACTION_PYTHON.is_file():
                venv.EnvBuilder(system_site_packages=True, with_pip=False).create(environment)

        def run_stage1(*arguments):
            status = system(shlex.join([str(EXTRACTION_PYTHON), '-u', *map(str, arguments)]))
            if status != 0:
                raise RuntimeError(f'Command stopped (exit {status}). See the output above.')

        run_stage1('run_layerwise_auroc.py', '--setup-only')
        """),
            md("""## 3. Authenticate and restore the pooled cache

        Enable notebook access to `HF_TOKEN_PLM` in Colab Secrets. The configured dataset
        must already exist and be private. No repository is created. A complete valid
        local cache skips downloads; a new runtime retrieves only the pooled arrays
        (about 516 MB before compression), not the roughly 63.5 GB residue maps."""),
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
        OUTPUT = Path('/content/nesg_esm_stage1_scratch/activations')
        OUTPUT.mkdir(parents=True, exist_ok=True)
        run_stage1('hf_storage.py', 'restore', '--repo', DESTINATION_REPO, '--folder', OUTPUT)
        assert (OUTPUT / 'completed.json').is_file(), 'The complete pooled extraction must be available in HF.'
        """),
            md("""## 4. Verify the development population

        IDs, original sequences and hashes are checked against the frozen inputs.
        Training/validation homology groups must be disjoint. Displayed counts are
        calculated from the local development manifest; no test records are loaded."""),
            code("""
        verify_inputs = (
            "import pandas as pd; from common import read_local; d=read_local(); "
            "assert d.split.value_counts().to_dict()=={'train':2361,'validation':691}; "
            "assert set(d.loc[d.split=='train','audit_cluster']).isdisjoint(set(d.loc[d.split=='validation','audit_cluster'])); "
            "print(pd.crosstab(d.split,d.label).rename(columns={0:'Insoluble',1:'Soluble'}).to_string()); "
            "print('Frozen sequence hashes and homology separation verified.')"
        )
        run_stage1('-c', verify_inputs)
        """),
            md("""## 5. Fit or reuse the 33 layer probes

        Each layer uses the same 1,280-dimensional mean-pooled representation and
        StandardScaler → L2 logistic regression pipeline. Scaling is learned inside
        each training fold. C is chosen by mean AUROC over the five frozen homology
        folds from `[0.0001, 0.001, 0.01, 0.1, 1, 10]`, then the final probe is refitted
        on all training proteins. Validation is never used to tune C or fit a model.

        The diagnostic composition baseline uses the 20 standard amino-acid fractions
        plus log length. X residues are retained without substitution. There is no
        class/family balancing. Validation uncertainty uses 2,000 homology-group
        bootstrap resamples shared across layers.

        Completed results are saved to private HF. Per-protein experimental labels
        remain in the separate local manifest and are omitted from published prediction
        tables. Fitted probes, coefficients, scores, plots and provenance are retained."""),
            code("""
        run_stage1(
            'run_layerwise_auroc.py',
            '--activations', OUTPUT,
            '--repo', DESTINATION_REPO,
            '--threads', '2',
        )
        latest = json.loads(Path('/content/nesg_layerwise_probes/latest.json').read_text())
        results = Path(latest['results'])
        print('Verified result directory:', results)
        """),
            md("""## 6. Training cross-validation versus fixed validation

        The train curve averages five held-out-fold AUROCs; its band is **±1 fold SD**.
        The validation curve uses the final model; its band is a **95% homology-group
        bootstrap CI**. These bands have different meanings. The same training folds
        select C, so train-CV scores are tuning estimates, not unbiased nested-CV results."""),
            code("""
        from IPython.display import Image, Markdown, display
        display(Image(filename=str(results / 'figures/train_validation_auroc.png')))
        """),
            md("""## 7. Training fit versus validation: overfitting diagnostic

        Here training AUROC is measured on the proteins used to fit the final probe.
        It is explicitly **in-sample** and can be optimistic. Use the previous figure
        for held-out-fold training estimates. A gap can reflect overfitting and
        differences between populations; it is not a causal diagnosis."""),
            code("""
        display(Image(filename=str(results / 'figures/train_fit_validation_auroc.png')))
        """),
            md("""## 8. Exact values, interpretation and saved outputs

        The largest validation AUROC is descriptive across 33 inspected layers;
        pointwise intervals do not correct for that selection. Layer-wise decoding
        does not establish direct versus indirect encoding or a causal computation.
        No SAE discovery, intervention or test evaluation is performed here."""),
            code("""
        import csv
        display(Markdown((results / 'TRAIN_VALIDATION_AUROC.md').read_text()))
        with (results / 'train_validation_auroc.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\\t'))
        lines = ['| Layer | Train CV | Train fit | Validation | Validation 95% CI |',
                 '|---:|---:|---:|---:|---:|']
        for row in rows:
            lines.append(
                f"| {row['layer']} | {float(row['cv_mean']):.4f} | {float(row['training_fit_auroc']):.4f} | "
                f"{float(row['auroc']):.4f} | {float(row['auroc_low']):.4f}–{float(row['auroc_high']):.4f} |"
            )
        display(Markdown('\\n'.join(lines)))
        url = f"https://huggingface.co/datasets/{latest['repository']}/tree/main/{latest['hf_prefix']}"
        display(Markdown(f'[Open private HF results]({url}) — download `results.zip` for PNG/SVG figures, exact tables and fitted probes.'))
        print('Plots are embedded above. No activations were recomputed.')
        """),
        ]
    )
    nb.metadata["colab"] = {"name": "ESM_Train_Validation_AUROC_Colab.ipynb"}
    return nb


if __name__ == "__main__":
    path = Path(__file__).resolve().parent / "ESM_Train_Validation_AUROC_Colab.ipynb"
    nb = build()
    for cell in nb.cells:
        if cell.cell_type == "code":
            compile(cell.source, "<notebook-cell>", "exec")
    nbformat.validate(nb)
    nbformat.write(nb, path)
    print(path)
