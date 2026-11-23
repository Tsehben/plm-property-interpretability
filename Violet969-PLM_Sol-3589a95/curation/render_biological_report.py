#!/usr/bin/env python3
"""Render the current auditable annotation snapshot; explicitly retain pending tools."""
import argparse,json
from pathlib import Path
import pandas as pd
from annotate_biology import ROOT,OUT,read,sha,write_json

def table(d):
    def fmt(x):
        if pd.isna(x):return '—'
        if isinstance(x,float):return f'{x:.4g}'
        return str(x)
    return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(fmt(v) for v in r)+' |' for r in d.itertuples(index=False,name=None))

def render(out):
    d=pd.read_csv(out/'protein_annotations.tsv',sep='\t',keep_default_na=False,na_values=['NA'])
    qc=json.loads((out/'validation_checks.json').read_text());summary=json.loads((out/'eda_summary.json').read_text())
    st=read(out/'annotation_status.tsv');coverage=st.groupby(['concept','status']).size().rename('n').reset_index()
    labels=d.groupby(['split','label']).size().rename('n').reset_index()
    fields=['full_length','full_charge_density_ph7p0','full_pi_estimate','full_hydropathy_mean','full_aggregation_fraction','full_disorder_fraction']
    descriptive=d[fields].agg(['mean','median','min','max']).T.reset_index(names='concept')
    support=read(out/'family_biochemical_support.tsv');support=support[support.n>=5]
    variation=read(out/'family_biochemical_variation.tsv')
    variation=variation.groupby('feature').within_sd_over_global_sd.median().reset_index()
    variance=read(out/'source_variance.tsv');variance=variance[variance.feature.isin(['full_hydropathy_mean','full_aggregation_fraction','full_disorder_fraction'])]
    corr=pd.read_csv(out/'concept_spearman.tsv',sep='\t',index_col=0)
    associations=read(out/'label_associations_source.tsv')
    examples=associations[(associations.center.isin(['JCSG','MCSG','SGPP']))&associations.feature.eq('full_disorder_fraction')][['split','center','n0','n1','standardized_difference']]
    tm_count=int(st[(st.concept=='transmembrane')&st.status.isin(['ok','zero_features'])].shape[0])
    tm_message=(f'Uniform DeepTMHMM2 annotation is complete for {tm_count:,} proteins.' if tm_count==13880 else
        f'**The authorized full DeepTMHMM2 CPU run is in progress. This merged snapshot contains {tm_count:,}/13,880 completed-shard predictions.** '
        'Check `transmembrane_progress.json` for live progress; the watcher refreshes this table and report on completion. The pilot estimated about 10–11 hours, with runtime depending on sequence lengths, startup and machine load.')
    export=json.loads((out/'server_upload_export.json').read_text()) if (out/'server_upload_export.json').exists() else {}
    report=f'''# UESolDS Phase 1: biological annotation and quality assessment

The annotation layer accounts for **all {len(d):,} proteins**, preserving the
heterogeneous A cohort, original labels/splits/sequences, stable IDs, hashes,
Pfam memberships and MMseqs2 groups. **The independent test was not opened.**
No solubility probe, SAE model, intervention or fine-tuning was performed.

**Current completion:** biochemical descriptors, the documented aggregation
reimplementation, learned disorder, reused Pfam/homology and uniform exact-sequence
taxonomy searches are complete. This is **not yet complete coverage of every
requested biological concept**: secondary structure, accessibility and SignalP
await the user's server outputs; the full TM job runs separately.

{tm_message}

## Cohort and QC

{table(labels)}

The proteins contain **{int(d.length.sum()):,} residues**; lengths range from
{int(d.length.min())} to {int(d.length.max())}. **{int(d.esm_truncated.sum())} proteins**
extend beyond the first 1,022 residues visible to the existing ESM activations.
Full and visible summaries are separate. No tag stripping, relabelling, balancing,
family exclusion, sequence replacement or new split was applied.

[Validation checks](biological_annotation/validation_checks.json) verify all
{qc['positional_records_checked']:,} positional records, unchanged input hashes,
{qc['domain_coordinates_checked']:,} retained Pfam domains with original scores,
and {qc['regional_intervals_checked']:,} interval coordinates. Scalar summaries
recomputed from saved positional data agree within {qc['scalar_position_max_error']:.3g}.
Detected homology groups remain disjoint across train/validation. Every sequence
is canonical in this frozen cohort; unknown/short-input handling is nevertheless
tested. No new biological value is filled from its solubility label.

## Annotation coverage and missingness

{table(coverage)}

`zero_features` is a successful prediction with no qualifying region; `no_hit`
is a successful Pfam search with no GA-passing domain. Taxonomic no-match,
ambiguous/engineered matches, pending computations and unavailable tools are
separate statuses. Missing numerical outputs remain `NA`, never zero.

Uniform Pfam coverage remains **{int(d.pfam_status.eq('assigned').sum()):,}/13,880**.
Taxonomy uses UniProt **2026_04** snapshots: exact archived sequence equality is
required. Species consensus covers **{int(d.taxonomy_status.eq('assigned_exact_species_consensus').sum()):,} proteins**;
there are **{int(d.taxonomy_status.eq('exact_match_ambiguous_taxonomy').sum()):,} ambiguous matches**,
**{int(d.taxonomy_status.eq('exact_match_engineered_taxon').sum())} engineered-taxon matches** and
**{int(d.taxonomy_status.eq('no_exact_match').sum()):,} no-exact-match records**.
A reference match does not establish an unmodified native construct or validate
the solubility label. One returned archive taxon lacked a taxonomy ID; it remains
ambiguous. No source/host/name inference filled this or any other gap.

Coverage tables stratify by [source × label × split](biological_annotation/coverage_source_label_split.tsv),
[Pfam × label × split](biological_annotation/coverage_family.tsv),
[length](biological_annotation/coverage_length.tsv) and
[terminal-motif flags](biological_annotation/coverage_construct.tsv).
Biochemical/aggregation/disorder computations have no input-dependent failures in
this cohort. Taxonomic missingness is reference coverage/ambiguity, not absence of
organism identity. Pending server tools cannot yet be assessed for source-dependent
prediction failures.

## Overall distributions and relationships

{table(descriptive)}

![Concept distributions](biological_annotation/concept_distributions.png)

Aggregation fraction is the fraction of residues in sequence-predicted APRs;
disorder fraction uses metapredict V3 scores ≥0.5. These are different constructs
and neither is an experimentally observed soluble fraction. Mean charge density
uses the stated pH and pKas, rather than the old KR−DE proxy.

Hydropathy and APR fraction have pooled Spearman correlation
**{corr.loc['full_hydropathy_mean','full_aggregation_fraction']:.3f}**;
APR and disorder fractions have correlation
**{corr.loc['full_aggregation_fraction','full_disorder_fraction']:.3f}**.
The aggregation scale shares hydrophobic/compositional information, so these
relationships are expected dependencies, not independent mechanistic validations.

![Concept relationships](biological_annotation/concept_relationships.png)

## Source and family variation

![Source-level biological variation](biological_annotation/source_concepts.png)

Sources occupy different descriptive biological profiles. The following fractions
of observed protein-level variance lie between source means, calculated separately
by original split. These are descriptive associations, not causal decompositions:

{table(variance)}

For family variation, {len(support):,} Pfam × split strata contain ≥5 proteins.
Their per-feature within-family SD relative to the pooled split SD has the
following median across strata:

{table(variation)}

![Within-family biochemical variation](biological_annotation/within_family_variation.png)

Thus family identity often narrows biochemical variation while leaving variation
for later conditional analysis. This is not a power calculation. A shared Pfam
domain does not hold the complete architecture, organism or experimental workflow
constant. Domain memberships overlap and repeated appearances are not independent
proteins. [Family support](biological_annotation/family_biochemical_support.tsv)
retains group/source counts and single-class families. Primary A is unchanged.

## Descriptive label associations

All prespecified descriptors and amino-acid fractions are retained irrespective
of association size. Separate tables give [pooled](biological_annotation/label_associations_pooled.tsv),
[within-source](biological_annotation/label_associations_source.tsv),
[within-family](biological_annotation/label_associations_family.tsv), and
[within-family/source](biological_annotation/label_associations_family_source.tsv)
contrasts, with ≥3 observations per label in each reported stratum and homology
group counts. These tables contain means and standardized differences, **not AUROC**.

For example, the label-1 minus label-0 disorder association varies by centre/split:

{table(examples)}

There are {summary['same_source_family_label_contrasts']} supported family/source/split
strata for these descriptive contrasts. This does not mean that all have matching
support in both splits; consult the earlier [family × source report](FAMILY_SOURCE_REPORT.md)
for controlled B feasibility. The previous 69-family focus is not a replacement
training cohort. Neither these associations nor the earlier provenance-sensitive
AUROCs establish direct versus indirect ESM encoding.

## Aggregation discrepancy and remaining limitations

The independent AGGRESCAN implementation reproduces the published SH3 example's
interior scores, APR counts and normalized means. The example's terminal values
use the reverse of the paper's stated terminal convention. We preserve the written
rule as primary and explicitly report the reversed-rule sensitivity.
**{int(d.aggregation_terminal_sensitivity_regions_changed.sum()):,} proteins** change
at least one APR boundary under that sensitivity. The mean absolute APR-fraction
change is **{(d.full_aggregation_fraction-d.aggregation_alternate_fraction).abs().mean():.5f}**;
the maximum is **{(d.full_aggregation_fraction-d.aggregation_alternate_fraction).abs().max():.4f}**.
This unresolved implementation ambiguity warrants exploratory use, particularly
for terminal localisation. **{int(d.aggregation_length_extrapolation.sum())} sequences**
exceed the old server's 2,000-residue limit and carry a separate extrapolation flag.

Metapredict is a validated accessible alternative, **not IUPred3**; its AlphaFold-
derived training targets introduce model dependence. DeepTMHMM2 and NetSurfP
depend directly on pretrained protein representations. Agreement with ESM features
will therefore not constitute independent experimental corroboration. Charge/pI
ignore folded microenvironments and modifications. Existing construct/tag flags
are observations, not confirmed engineering or truncation annotations. Pfam clan
competition is retained; candidate architectures may contain overlapping alternatives.
UESolDS label-provenance limitations remain unchanged.

## Readiness and next action

- **Ready as reproducible sequence concepts:** length, composition, charge estimates,
  pI and hydropathy, using ESM-visible summaries when relating to cached activations.
- **Ready with explicit annotation uncertainty:** Pfam memberships/architectures and
  homology groups. Preserve overlap/no-hit indicators and the prior controlled subsets.
- **Exploratory:** aggregation (including terminal sensitivity), learned disorder,
  exact archive species-consensus taxonomy and model-derived TM predictions once complete.
  They must not be described as experimental biology.
- **Pending:** NetSurfP structure/accessibility and SignalP results; no placeholder
  values or generic hydrophobicity substitutions are supplied. Primary A keeps all
  records while these annotations are missing.

Continue using heterogeneous all-source A with source-stratified and within-family
analyses. Annotation does not solve uncertain experimental labels or the sparse,
overlap-sensitive same-source B support established earlier. Do not select a cohort
or explanatory feature set based on solubility prediction performance.

The server packet is at `{export.get('archive','not exported')}`. It contains
20 NetSurfP batches covering all 13,880 records; SignalP-ready batches cover
3,014 Eukarya and 3,182 Bacteria/Archaea records, with 7,684 unresolved settings
retained separately. Follow its README and retain raw positional outputs/job metadata.
Full sequences and neutral upload IDs have reproducible hash mappings. Files are
outside the checkout to avoid temporary FASTAs following branch switches.

## Files

- [Methods, dependencies and reproduction](BIOLOGICAL_ANNOTATION_METHODS.md)
- [Unified protein table](biological_annotation/protein_annotations.tsv)
- [Column dictionary](biological_annotation/feature_dictionary.tsv)
- [Per-record annotation status](biological_annotation/annotation_status.tsv)
- [Regional coordinates](biological_annotation/regions.tsv) and [Pfam domains](biological_annotation/pfam_domains.tsv)
- `aggregation_profiles.npz`, `sequence_profiles.npz`, `disorder_shards/`: full positional data
- `taxonomy_responses/`, `taxonomy.tsv`: exact-match reference evidence
- `transmembrane_pilot_results/`, `transmembrane_shards/`: raw predictor evidence
- [Input hashes](biological_annotation/input_hashes.json), [QC](biological_annotation/validation_checks.json), [output hashes](biological_annotation/output_hashes.json)
'''
    target=ROOT/'BIOLOGICAL_ANNOTATION_REPORT.md' if out.resolve()==OUT.resolve() else out/'BIOLOGICAL_ANNOTATION_REPORT.md'
    target.write_text(report)
    write_json(out/'output_hashes.json',{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and
        p.name not in ['output_hashes.json','transmembrane_progress.json','transmembrane_refresh.log'] and
        'transmembrane_shards' not in p.parts})  # each completed TM shard carries its own immutable hashes
    print(str(target))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=OUT);a=p.parse_args();render(a.out)
