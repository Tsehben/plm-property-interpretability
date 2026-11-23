"""Create the separate evaluation-only notebook; never overwrite the running notebook."""
from pathlib import Path
import nbformat
from build_notebooks import md,code,notebook
from build_property_notebook import build as uncontrolled


def build():
    setup=uncontrolled().cells[1:7]
    setup[1].source=setup[1].source.replace("/content/nesg_property_probes","/content/nesg_within_eval_code")
    setup[1].source=setup[1].source.replace("'property_inputs.py', 'property_probes.py',", "'within_property_evaluation.py', 'run_within_property.py', 'property_inputs.py', 'property_probes.py',")
    setup[3].source=setup[3].source.replace("EXTRACTION_PYTHON","EVALUATION_PYTHON")
    setup[4]=md("""## 3. Authenticate for saved probe results

Use the existing private HF repository and `HF_TOKEN_PLM` Colab Secret. The next stage downloads
saved probe checkpoints only when they are not already available in this runtime. The checkpoint
archives contain predictions and fitted-probe files; only prediction tables are extracted for use.
No ESM weights, pooled activations or residue maps are downloaded. No repository is created.""")
    setup[5].source=setup[5].source.replace("OUTPUT = Path('/content/nesg_esm_stage1_scratch/activations')\nOUTPUT.mkdir(parents=True, exist_ok=True)\n","")
    setup[5].source=setup[5].source.replace("print('Configured private HF destination. Cached results will be checked before downloading activations.')", "print('Configured private HF destination for saved probe predictions only.')")
    nb=notebook([
        md("""# NESG: within-property evaluation of existing probes

**Let the uncontrolled notebook finish and inspect its results first. Then run this notebook.**
Keep all **33 layers**: the expensive probe fitting has already been done. This notebook recalculates
metrics from saved validation predictions inside biological subsets. It neither fits new probes nor
runs protein sequences through ESM. A CPU runtime is sufficient.

| Question | Evaluation |
|---|---|
| Is soluble recovery still distinguishable within a hydrophobicity range? | Same global solubility probe, evaluated separately in training-defined ranges |
| What about charge, length, aggregation or APR burden? | Restrict one property at a time using the same fixed rule |
| Is solubility distinguishable within a family? | Same global solubility probe within each observed Pfam family; support shown explicitly |
| Are other properties still decodable in these subsets? | Reuse continuous-property predictions within ranges/families and family-membership predictions within numeric ranges |

There is **no prediction of the defining family's membership inside its own family**, where that target
is constant. No joint matching, residualization, retraining on subsets, SAE analysis, interventions or
reserved-test evaluation is performed. Broad range restriction is an approximate restriction, not
proof that a confounder has been removed. The canonical heterogeneous dataset is unchanged.

**Status:** input/support and engineering checks are performed locally. Actual all-layer conditional
performance appears only after this notebook runs against your completed private HF results.
The current uncontrolled notebook and its checkpoints are preserved."""),*setup,
        md("""## 4. Define restrictions without using probe performance

For each continuous descriptor, the four cut points are training-only 20th/40th/60th/80th percentiles.
Tied boundaries are collapsed; lower bounds are inclusive, upper bounds exclusive, and the outer
intervals extend to infinity. Missing descriptors are never imputed. These are **one-property-at-a-time**
restrictions; other biological properties remain heterogeneous.

Families come from the existing all-GA Pfam assignments observed in training. No ≥3/≥5 training-solubility
threshold is imposed because no new within-family models are being fitted. Support tables include every
training-observed family. Point estimates require at least two validation proteins and, for AUROC/AP,
both labels. Singleton/single-class families remain explicitly unsupported. Proteins may appear in
several distinct family analyses but are counted only once within each family.

The noncompeting-assignment sensitivity removes flagged memberships only for that family analysis;
it does not change the primary cohort. Both variants retain their own denominators."""),
        code("""
        run_stage1('-c', r'''
from property_inputs import load_inputs
from within_property_evaluation import plan
d,p,m,c=load_inputs()
targets,jobs,support,edges,families=plan(d,p,m,c)
print('Frozen split sizes:', d.split.value_counts().to_dict())
print('Training-only cut points:')
for name, cuts in edges.items(): print(name, cuts)
print('\\nSolubility support (not performance):')
print(support[support.target.eq('solubility')].groupby(['kind','variant','status']).size().to_string())
print('\\nNo prediction files loaded; no fitting; reserved test untouched.')
''')
        """),
        md("""## 5. Evaluate all layers using saved predictions

The source must contain all 33 completed uncontrolled layers plus composition and its final summary.
If that run is still underway, this cell stops with a message to let it finish; it does not interrupt
or refit it. Leave `SOURCE_RUN` empty to use the first compatible completed run in deterministic prefix
order. You may set the exact `esm_stage1/uncontrolled_property_probing_v1/...` prefix to choose an
existing run explicitly. Selection never uses performance.

ID/hash/split alignment and target-specific coverage are checked against frozen local manifests.
The same global predictions are reused for every subset. Composition comparisons use exactly the
same proteins within a subset. No train-fit performance is substituted for validation performance.

Solubility uses 2,000 shared homology-group bootstrap draws. Intervals require ≥3 groups per label
and ≥90% valid draws. These are minimum reporting safeguards, not confirmatory-power thresholds.
Secondary targets have point estimates only, keeping this evaluation lighter than fitting new models.
Completed local metric calculations resume by layer; final results are saved to private HF. Losing
the runtime before final upload can require recomputing metrics, but never refitting the probes."""),
        code("""
        SOURCE_RUN = ''
        arguments = ['run_within_property.py', '--repo', DESTINATION_REPO]
        if SOURCE_RUN:
            arguments += ['--source-run', SOURCE_RUN]
        run_stage1(*arguments)
        latest = json.loads(Path('/content/nesg_within_property/latest.json').read_text())
        results = Path(latest['results'])
        from IPython.display import Image, Markdown, HTML, display
        print('Verified evaluation results:', results)
        """),
        md("""## 6. Solubility within numeric property ranges"""),
        code("display(Image(filename=str(results / 'figures/solubility_within_ranges.png')))"),
        md("""**Interpretation:** each colored curve uses the same saved global probe on a different
validation subset. Read its soluble/insoluble counts and interval before interpreting its level or
shape. The dashed overall curve uses a different population; its difference from a range-specific
curve is descriptive. A drop can reflect reduced biological variation, changed difficulty, changed
label mix or sampling uncertainty. It is not a causal estimate of reliance on that property.

All 33 layers stay in the analysis. Choosing 4/16/32 later for a compact illustration does not require
discarding other layers. Do not choose layers solely because their restricted AUROC drops the most."""),
        md("""## 7. Within-family solubility and competing-assignment sensitivity"""),
        code("""
        display(Image(filename=str(results / 'figures/solubility_within_families_all_ga.png')))
        display(Image(filename=str(results / 'figures/solubility_within_families_noncompeting.png')))
        """),
        md("""**Interpretation:** family rows are separate, overlapping evaluations, not independent
replicates. One soluble and one insoluble protein can produce AUROC 0 or 1 with almost no inferential
support. The labels show protein and independent-group counts precisely to prevent treating such
rows as strong evidence. Competing-hit removal changes the evaluation population and can leave very
few contrasts. Shared Pfam membership does not establish equal architecture, organism or assay."""),
        code("""
        import csv, math, html
        with (results / 'stratum_target_support.tsv').open() as handle:
            support_rows = list(csv.DictReader(handle, delimiter='\\t'))
        chosen = [r for r in support_rows if r['target']=='solubility' and r['kind']=='family' and r['status']=='estimable']
        columns = ['control','category','variant','train_n','validation_n','validation_positive','positive_groups','negative_groups','ci_support_eligible']
        lines = ['| ' + ' | '.join(columns) + ' |', '| ' + ' | '.join(['---']*len(columns)) + ' |']
        lines += ['| ' + ' | '.join(r.get(c,'') for c in columns) + ' |' for r in chosen]
        display(Markdown('\\n'.join(lines)))
        """),
        md("""## 8. Other properties within ranges

These panels use continuous-target R², not AUROC. Full MAE/RMSE results and within-family continuous
results are in the exact-result viewer below. Family-membership AUROC/AP inside numeric ranges is
also available there. The subset definitions are identical across layers for each target."""),
        code("""
        for name in ['hydrophobicity','charge_proxy','log_length','aggregation_Na4vSS','apr_fraction']:
            display(Image(filename=str(results / f'figures/secondary_{name}.png')))
        """),
        md("""**Interpretation:** restricting the property being predicted reduces target variance,
which can sharply reduce R² even if absolute error changes little. Inspect target SD and MAE together.
Constant targets have undefined R², not a meaningful failure to decode. Hydropathy, charge and length
are deterministic descriptors; AGGRESCAN outputs are predictions rather than aggregation measurements.
These secondary point estimates have no confidence bands and are exploratory."""),
        md("""## 9. Exact results for any target/subset

Change the two display strings below and rerun **only this cell** to inspect a different subset.
Examples: `range:hydrophobicity:1`, `range:aggregation_Na4vSS:3`, `family:PF04055:all_ga`,
`family:PF04055:noncompeting`, or `overall`. `TARGET` can be `solubility`, one of the five continuous
target names, or an eligible Pfam accession for numeric-range evaluation. This changes the display only."""),
        code("""
        TARGET = 'solubility'
        STRATUM = 'range:hydrophobicity:3'
        with (results / 'within_property_metrics.tsv').open() as handle:
            rows = [r for r in csv.DictReader(handle, delimiter='\\t') if r['target']==TARGET and r['key']==STRATUM]
        if not rows:
            print('No evaluated task for this combination. Check stratum_target_support.tsv; family-membership targets are evaluated only in numeric ranges.')
        else:
            first = rows[0]
            display(Markdown(f"**{TARGET} / {STRATUM}** — validation N={first['validation_n']}; groups={first['validation_groups']}; target SD={first['validation_target_sd']}; status={first['status']}"))
            columns = ['layer','metric','value','low','high','composition','delta_vs_composition','delta_low','delta_high','bootstrap_valid']
            def formatted(value):
                try:
                    v=float(value)
                    return f'{v:.4f}' if math.isfinite(v) else 'undefined'
                except (ValueError,TypeError):
                    return value or 'undefined'
            body='<tr>'+''.join('<th>'+html.escape(c)+'</th>' for c in columns)+'</tr>'
            for row in rows:
                body+='<tr>'+''.join('<td>'+html.escape(formatted(row[c]))+'</td>' for c in columns)+'</tr>'
            display(HTML('<div style="max-height:600px;overflow:auto"><table>'+body+'</table></div>'))
        """),
        md("""## 10. Saved report and boundaries of the analysis"""),
        code("""
        display(Markdown((results / 'WITHIN_PROPERTY_RESULTS.md').read_text()))
        url=f"https://huggingface.co/datasets/{latest['repository']}/tree/main/{latest['hf_prefix']}"
        display(Markdown(f'[Private HF evaluation results]({url}) — `results.zip` contains the complete tables, protocol, source hashes and PNG/SVG figures.'))
        """),
        md("""Keep the complete layer-wise curves as the primary descriptive result. Small differences
or the largest observed drop across many restrictions need uncertainty and multiplicity caution.
This notebook does not establish a direct solubility representation or a causal computation.
Mechanistic hypotheses and SAE interventions remain separate subsequent work."""),
    ])
    nb.metadata["colab"]={"name":"ESM_Within_Property_Evaluation_Colab.ipynb"}
    return nb


if __name__=="__main__":
    nb=build()
    for cell in nb.cells:
        if cell.cell_type=="code": compile(cell.source,"<cell>","exec")
    nbformat.validate(nb)
    p=Path(__file__).resolve().parent/"ESM_Within_Property_Evaluation_Colab.ipynb"
    nbformat.write(nb,p)
    print(p)
