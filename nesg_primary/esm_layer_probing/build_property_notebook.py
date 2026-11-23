"""Build a separate uncontrolled-property notebook; preserve existing notebooks."""
from pathlib import Path
import nbformat
from build_notebooks import md, code, notebook
from build_auroc_notebook import build as old_build


def build():
    setup = old_build().cells[1:7]
    setup[1].source = setup[1].source.replace("work = Path('/content/nesg_esm_stage1')", "work = Path('/content/nesg_property_probes')")
    setup[1].source = setup[1].source.replace("'data/input_integrity.json', 'data/local_labels.tsv',", """'data/input_integrity.json', 'data/local_labels.tsv',
    'property_inputs.py', 'property_probes.py', 'property_report.py', 'run_property_probes.py',
    'data/pfam_memberships.tsv', 'data/uncontrolled/properties.tsv',
    'data/uncontrolled/family_catalog.tsv', 'data/uncontrolled/protocol.json',
    'data/uncontrolled/coverage.tsv', 'data/uncontrolled/integrity.json',""")
    setup[4] = md("""## 3. Authenticate securely

Enable notebook access to `HF_TOKEN_PLM` in Colab Secrets. The configured dataset must already
exist and be private. No repository is created. The runner first checks saved results using small
metadata files; it downloads pooled arrays only if an unfinished layer needs fitting. Completed
runs reopen without downloading activations or refitting. Full residue maps are never downloaded.""")
    setup[5].source = setup[5].source.replace(
        "run_stage1('hf_storage.py', 'restore', '--repo', DESTINATION_REPO, '--folder', OUTPUT)\nassert (OUTPUT / 'completed.json').is_file(), 'The complete pooled extraction must be available in HF.'",
        "print('Configured private HF destination. Cached results will be checked before downloading activations.')")
    cells = [md("""# NESG: uncontrolled property decoding across all 33 ESM-2 layers

**Run this notebook first. Review its performance plots before designing controlled comparisons.**
Choose **Run all** in Colab; a CPU runtime is sufficient. No ESM extraction, SAE analysis,
restriction/matching, or reserved-test evaluation occurs here. The heterogeneous cohort
and original family distributions remain unchanged. All 33 ESM layers are used.

| Target | Meaning | Probe | Evaluation |
|---|---|---|---|
| Solubility | NESG experimental soluble recovery at high expression | L2 logistic; reuse completed fits | AUROC, AP |
| Hydrophobicity | Mean Kyte–Doolittle sequence descriptor | Ridge | R², MAE, RMSE |
| Charge proxy | (K+R−D−E)/length, not pH-dependent charge | Ridge | R², MAE, RMSE |
| Length | Natural log original sequence length | Ridge | R², MAE, RMSE |
| Aggregation propensity | AGGRESCAN Na4vSS prediction | Ridge | R², MAE, RMSE |
| APR burden | Predicted fraction of residues in aggregation-prone regions | Ridge | R², MAE, RMSE |
| Pfam membership | One presence/absence target per training-eligible family | Separate L2 logistic probes | Per-family and macro AUROC/AP |

**Persistence:** pooled activations and result checkpoints live in private HF. Drive provides code
and frozen inputs only. Complete matching runs skip fitting; partial runs resume completed layers
from HF and completed targets within the same Colab runtime. No model weights or full maps are downloaded.

**Execution status:** actual composition fits and numerical/integrity tests were run locally.
The new all-layer property results require your authenticated Colab run; there are no placeholder curves.
Existing solubility curves are reused when a compatible completed bundle is available."""), *setup,
        md("""## 4. Verify annotations and target eligibility

The frozen manifests contain **2,361 training and 691 validation proteins**. The three sequences
containing X remain in the cohort. For hydropathy, charge proxy and AGGRESCAN targets, their undefined
values are omitted only from that target's fits/evaluation: **2,359 training and 690 validation**.
Zero APR burden is a valid prediction and is retained.

Family targets require ≥5 positive and ≥5 negative **training homology groups**, with both membership
classes in every held-out fold and fitting complement. This yields **21 families**; it selects target
tasks, not proteins. Every protein remains an example for every family task. A protein can be positive
for several Pfam families. This is family-membership prediction, not within-family solubility evaluation.
Absence of a GA-passing hit is an annotation-defined negative, not proof of biological absence.
Of these 21 training-eligible targets, **16 have positives in validation**; the other five remain
in the output with undefined discrimination. Some of the 16 have only one validation positive."""),
        code("""
        run_stage1('-c', r'''
from property_inputs import load_inputs
from pathlib import Path
import pandas as pd
d,p,m,c=load_inputs()
assert d.split.value_counts().to_dict()=={'train':2361,'validation':691}
assert set(d.loc[d.split=='train','audit_cluster']).isdisjoint(d.loc[d.split=='validation','audit_cluster'])
print(pd.read_csv('data/uncontrolled/coverage.tsv',sep='\\t').to_string(index=False))
print('\\nFamily targets selected using training only:')
print(c.loc[c.eligible,['pfam','name','train_positive','train_negative','train_positive_groups']].to_string(index=False))
print('\\nFrozen sequence hashes, annotation hashes, splits and homology separation passed.')
''')
        """),
        md("""## 5. Fit or resume uncontrolled property probes

Standardization is fitted separately inside each training fold. Binary probes select C from
`[0.0001, 0.001, 0.01, 0.1, 1, 10]` by mean fold AUROC. Continuous probes select ridge alpha from
`[0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000]` by minimum mean fold MSE. Targets remain in
their original units (length is explicitly log-transformed). Each selected model is refitted on all
observed training examples and evaluated on the fixed validation examples for that target.

There is no class balancing or validation-driven tuning. Each target has a composition + log-length
baseline and a training-mean/prevalence constant reference. The same training folds tune and report
CV performance; these are tuning estimates, **not nested-CV generalization estimates**.

Completed solubility fits are matched by input/protocol/code/activation hashes and reused. If none
match, the runner announces that it will fit them under the unchanged protocol. Every completed
layer is uploaded atomically with hashes to the existing private HF dataset. No observed per-protein
targets or sequences are included in uploaded prediction tables. IDs/hashes permit local annotation joins.

This fits many more models than the solubility-only notebook. Progress prints each completed target;
if interrupted, rerun this cell. A deleted runtime loses at most the unfinished layer's fits."""),
        code("""
        run_stage1('run_property_probes.py', '--activations', OUTPUT,
                   '--repo', DESTINATION_REPO, '--threads', '2')
        latest = json.loads(Path('/content/nesg_uncontrolled_probes/latest.json').read_text())
        results = Path(latest['results'])
        from IPython.display import Image, Markdown, display
        print('Verified results:', results)
        """),
        md("""## 6. Experimental soluble-recovery decoding

The train curve is mean held-out-fold AUROC (±1 SD); validation has a pointwise 95% homology-group
bootstrap CI. The two bands have different meanings. Existing validated fits are reused, so this
panel connects the new property analysis to the previous solubility result."""),
        code("display(Image(filename=str(results / 'figures/solubility_auroc.png')))"),
        md("""**Interpretation:** above-chance discrimination shows linearly accessible association with
NESG soluble recovery. It does not tell us whether the model uses family, composition, aggregation,
or other features, and does not establish a causal mechanism."""),
        md("""## 7. Continuous biological descriptors and predicted aggregation"""),
        code("display(Image(filename=str(results / 'figures/continuous_r2.png')))"),
        md("""**Interpretation:** R² = 1 indicates perfect prediction; R² = 0 matches the validation-mean
reference, and negative values are possible. The separate training-mean constant baseline is saved
in the tables and need not have R² exactly zero on validation. Hydrophobicity and charge proxy are
linear functions of composition, while log length is directly supplied to the baseline. Near-perfect
composition performance is therefore expected. ESM performance measures retention/linear accessibility
of these descriptors, not additional experimental knowledge. AGGRESCAN outputs are predictions."""),
        code("display(Image(filename=str(results / 'figures/continuous_mae.png')))"),
        md("""**Interpretation:** MAE is in each target's units and lower is better. Compare layers within
a target; do not rank different properties by their raw MAEs or compare R² numerically with AUROC.
Only paired differences measured on the same target population support ESM-versus-composition comparisons."""),
        md("""## 8. Pfam family-membership decoding"""),
        code("""
        display(Image(filename=str(results / 'figures/family_macro.png')))
        display(Image(filename=str(results / 'figures/family_heatmap.png')))
        """),
        md("""**Interpretation:** each family is a separate membership target, with all heterogeneous
proteins retained. Macro averages give each evaluable family equal weight. AP must be read against
family prevalence. A family without both membership classes in validation has undefined discrimination
and is reported explicitly. The same evaluable family set is used across layers and valid bootstrap
draws. Sparse positive support can make intervals wide or reduce valid bootstrap draws.
Macro intervals are omitted if fewer than 90% of draws support every included family; exact valid
draw counts are saved in `family_macro.tsv`. A missing band is not evidence of zero uncertainty.

These are annotation-presence probes. All GA-passing memberships, including competing overlapping
assignments, are retained. This notebook does not resolve architecture ambiguity or establish
within-family solubility contrasts. Predicting a family's membership *inside a subset where every
protein has that family* would be a constant-target task and is not a useful controlled probe."""),
        md("""## 9. Exact results and the decision before controlled analyses"""),
        code("""
        import csv
        def show_tsv(path, columns):
            with path.open() as handle:
                rows = list(csv.DictReader(handle, delimiter='\\t'))
            text = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---']*len(columns)) + ' |']
            for row in rows:
                text.append('| ' + ' | '.join(str(row.get(c,'')) for c in columns) + ' |')
            display(Markdown('\\n'.join(text)))
        show_tsv(results / 'family_support.tsv', ['pfam','name','train_positive','validation_positive','validation_negative','validation_estimable'])
        from IPython.display import HTML
        import html, math
        with (results / 'metrics.tsv').open() as handle:
            all_metrics = list(csv.DictReader(handle, delimiter='\\t'))
        for target in dict.fromkeys(row['target'] for row in all_metrics):
            rows = [r for r in all_metrics if r['target'] == target and r['metric'] in ('r2','auroc')]
            def number(value):
                try:
                    x = float(value)
                    return f'{x:.4f}' if math.isfinite(x) else 'undefined'
                except (ValueError, TypeError):
                    return 'undefined'
            body = '<tr><th>Layer (0 = composition)</th><th>Metric</th><th>Train CV</th><th>Validation</th><th>95% CI</th><th>Valid bootstrap draws</th></tr>'
            for row in rows:
                cells = [row['layer'], row['metric'], number(row['cv_mean']), number(row['validation']),
                         number(row['low']) + ' – ' + number(row['high']), row['bootstrap_valid']]
                body += '<tr>' + ''.join('<td>' + html.escape(x) + '</td>' for x in cells) + '</tr>'
            display(HTML('<details><summary>' + html.escape(target) + ': exact results for all layers</summary><table>' + body + '</table></details>'))
        display(Markdown((results / 'UNCONTROLLED_RESULTS.md').read_text()))
        url = f"https://huggingface.co/datasets/{latest['repository']}/tree/main/{latest['hf_prefix']}"
        display(Markdown(f'[Private HF checkpoints and figures]({url}). Summary tables/plots are in `summary/results.zip`; each `layer_XX/results.zip` contains fitted parameters, CV results and predictions. Layer 00 is composition.'))
        """),
        md("""**Stop here and inspect the results.** Next we can prespecify controlled comparisons using
the same globally trained solubility/property probes: for example within supported families or within
training-defined hydrophobicity ranges, provided target variation and independent validation support
remain. Restricted-population AUROC changes alone are not causal estimates of confounding and can
reflect changed difficulty, target range and sample size. Family membership is constant inside one
family; it should not be evaluated as a discrimination task there. Do not choose restrictions because
they yield the largest observed AUROC drop. No such analysis is run by this notebook."""),
    ]
    nb = notebook(cells)
    nb.metadata["colab"] = {"name": "ESM_Uncontrolled_Property_Probes_Colab.ipynb"}
    return nb


if __name__ == "__main__":
    nb = build()
    for cell in nb.cells:
        if cell.cell_type == "code":
            compile(cell.source,"<notebook-cell>","exec")
    nbformat.validate(nb)
    path = Path(__file__).resolve().parent/"ESM_Uncontrolled_Property_Probes_Colab.ipynb"
    nbformat.write(nb,path)
    print(path)
