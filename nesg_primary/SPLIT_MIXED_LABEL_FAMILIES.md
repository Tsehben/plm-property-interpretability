# Split-specific mixed-label Pfam families

**Held-out-label inspection disclosure.** The user explicitly requested this analysis of reserved-test labels. Producing and viewing the test-label-by-family tables constitutes inspection of held-out labels. The test can no longer be described as uninspected with respect to family–label composition. This is not model-performance evaluation, and these results must not be used to revise the cohort, split assignment, model selection or experimental hypotheses. No such changes were made. The UESolDS independent test was not accessed.

The original `results/test_lock.json` is preserved as a historical freeze record; its unchanged `LOCKED` field does not negate this subsequent authorized inspection. Access restrictions for ordinary notebook reruns remain in place. This standalone inspection is recorded in `TEST_LABEL_EXPOSURE.md` and its machine-readable inspection event.

## Definitions and validation

- Each split is analyzed independently. A qualifying family has ≥1 soluble and ≥1 insoluble protein **in that split**; no cross-split support rule, ≥3/≥5 threshold or model result is used. Training is evaluated independently by the same definition.
- Soluble is label 1; insoluble is label 0. Protein counts deduplicate stable record ID and Pfam accession. Each split has unique sequence hashes, so these also equal unique sequences. Repeat domains/assignments do not inflate counts.
- Family names are the Pfam 37.0 biological descriptions. Full CSV/TSV files also contain the short Pfam name. These are domain families, not necessarily identical full-protein architectures.
- Homology counts use the existing audit MMseqs2 component IDs. A group containing both labels counts once in each class; the CSV/TSV tables additionally report shared groups. Families, clans and homology groups remain distinct.
- *Competing overlap* means a hit for this family overlaps a different-family hit from the same nonempty Pfam clan on the same protein (1-based inclusive aligned coordinates). Existing flags were checked directly against coordinates.
- The primary tables retain **all GA-passing assignments**. The separate sensitivity removes both competing hits, then deduplicates memberships and independently reassesses mixed status. An unflagged repeat of the same family can preserve membership. No proteins or labels are removed from the cohort.
- *Affected proteins* counts unique proteins in that table with an original competing hit for the reported family. Sensitivity-table flags retain the original evidence among surviving memberships; they do not imply that flagged hits were retained.
- Tables are complete, sorted by total unique proteins descending, then accession. Summing rows can count multidomain proteins more than once; split-summary protein totals use unions.
- Input SHA256 matches the frozen manifests and pinned domain artifact. Split sizes/labels, unique IDs/hashes, coordinate bounds, sequence lengths, both GA thresholds, overlap flags and absence of group/hash overlap between every pair of splits were asserted. Input hashes are unchanged after execution.

## Counts

| Split | Split proteins (soluble / insoluble) | All-GA mixed families | Mixed families after overlap removal |
| --- | ---: | ---: | ---: |
| train | 2361 (1645 / 716) | 237 | 150 |
| validation | 691 (491 / 200) | 81 | 38 |
| reserved_test | 710 (478 / 232) | 100 | 48 |

## Primary: all GA-passing assignments

### train: 237 qualifying families

[Complete CSV](results/split_mixed_label_families/train_all_GA_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/train_all_GA_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF13450 | NAD(P)-binding Rossmann-like domain | 28 | 41 | 69 | 11 | 26 | Yes | 67 |
| PF01266 | FAD dependent oxidoreductase | 18 | 44 | 62 | 6 | 21 | Yes | 47 |
| PF00455 | DeoR C terminal sensor domain | 31 | 26 | 57 | 1 | 1 | No | 0 |
| PF08220 | DeoR-like helix-turn-helix domain | 32 | 25 | 57 | 2 | 1 | Yes | 41 |
| PF04055 | Radical SAM superfamily | 31 | 21 | 52 | 4 | 11 | No | 0 |
| PF00575 | S1 RNA binding domain | 48 | 2 | 50 | 6 | 1 | Yes | 2 |
| PF00535 | Glycosyl transferase family 2 | 37 | 11 | 48 | 24 | 5 | Yes | 33 |
| PF01494 | FAD binding domain | 19 | 26 | 45 | 6 | 13 | Yes | 34 |
| PF01656 | CobQ/CobB/MinD/ParA nucleotide binding domain | 31 | 11 | 42 | 3 | 5 | Yes | 39 |
| PF00890 | FAD binding domain | 11 | 29 | 40 | 8 | 19 | Yes | 39 |
| PF08279 | HTH domain | 24 | 16 | 40 | 2 | 1 | Yes | 40 |
| PF00919 | Uncharacterized protein family UPF0004 | 25 | 11 | 36 | 1 | 1 | No | 0 |
| PF13649 | Methyltransferase domain | 18 | 17 | 35 | 8 | 11 | Yes | 35 |
| PF08241 | Methyltransferase domain | 18 | 16 | 34 | 8 | 10 | Yes | 34 |
| PF07992 | Pyridine nucleotide-disulphide oxidoreductase | 15 | 17 | 32 | 6 | 12 | Yes | 26 |
| PF13614 | AAA domain | 29 | 1 | 30 | 2 | 1 | Yes | 30 |
| PF13641 | Glycosyltransferase like family 2 | 21 | 9 | 30 | 15 | 4 | Yes | 30 |
| PF03061 | Thioesterase superfamily | 21 | 8 | 29 | 7 | 2 | Yes | 23 |
| PF09140 | ATPase MipZ | 28 | 1 | 29 | 2 | 1 | Yes | 29 |
| PF10609 | NUBPL iron-transfer P-loop NTPase | 27 | 2 | 29 | 2 | 2 | Yes | 29 |
| PF07282 | Putative transposase DNA-binding domain | 22 | 6 | 28 | 10 | 5 | No | 0 |
| PF01636 | Phosphotransferase enzyme family | 14 | 13 | 27 | 4 | 6 | Yes | 13 |
| PF01966 | HD domain | 23 | 4 | 27 | 3 | 2 | Yes | 3 |
| PF08242 | Methyltransferase domain | 13 | 13 | 26 | 8 | 8 | Yes | 26 |
| PF13489 | Methyltransferase domain | 12 | 14 | 26 | 6 | 8 | Yes | 26 |
| PF01926 | 50S ribosome-binding GTPase | 20 | 4 | 24 | 3 | 1 | Yes | 16 |
| PF01938 | TRAM domain | 17 | 7 | 24 | 2 | 1 | Yes | 2 |
| PF01588 | Putative tRNA binding domain | 20 | 2 | 22 | 2 | 1 | No | 0 |
| PF12831 | FAD dependent oxidoreductase | 5 | 16 | 21 | 4 | 12 | Yes | 19 |
| PF13847 | Methyltransferase domain | 12 | 9 | 21 | 7 | 7 | Yes | 21 |
| PF02374 | Anion-transporting ATPase | 16 | 4 | 20 | 3 | 2 | Yes | 16 |
| PF07728 | AAA domain (dynein-related subfamily) | 5 | 14 | 19 | 4 | 8 | Yes | 6 |
| PF10111 | Glycosyltransferase like family 2 | 13 | 6 | 19 | 11 | 3 | Yes | 19 |
| PF01368 | DHH family | 16 | 2 | 18 | 1 | 1 | No | 0 |
| PF01385 | Probable transposase | 17 | 1 | 18 | 5 | 1 | No | 0 |
| PF13167 | GTP-binding GTPase N-terminal | 14 | 4 | 18 | 1 | 1 | No | 0 |
| PF16360 | GTP-binding GTPase Middle Region | 14 | 4 | 18 | 1 | 1 | No | 0 |
| PF00004 | ATPase family associated with various cellular activities (AAA) | 7 | 10 | 17 | 6 | 6 | Yes | 16 |
| PF02272 | DHHA1 domain | 14 | 2 | 16 | 1 | 1 | No | 0 |
| PF03966 | Trm112p-like protein | 15 | 1 | 16 | 2 | 1 | No | 0 |
| PF04390 | Lipopolysaccharide-assembly | 15 | 1 | 16 | 2 | 1 | No | 0 |
| PF00696 | Amino acid kinase family | 11 | 4 | 15 | 2 | 3 | No | 0 |
| PF01593 | Flavin containing amine oxidoreductase | 4 | 11 | 15 | 4 | 7 | Yes | 14 |
| PF16901 | C-terminal domain of alpha-glycerophosphate oxidase | 10 | 5 | 15 | 1 | 1 | No | 0 |
| PF02666 | Phosphatidylserine decarboxylase | 9 | 5 | 14 | 2 | 3 | No | 0 |
| PF04012 | PspA/IM30 family | 7 | 7 | 14 | 1 | 1 | No | 0 |
| PF05175 | Methyltransferase small domain | 8 | 6 | 14 | 4 | 5 | Yes | 14 |
| PF13353 | 4Fe-4S single cluster domain | 10 | 4 | 14 | 2 | 2 | No | 0 |
| PF13622 | Acyl-CoA thioesterase N-terminal domain | 8 | 6 | 14 | 2 | 1 | Yes | 14 |
| PF13738 | Pyridine nucleotide-disulphide oxidoreductase | 8 | 6 | 14 | 1 | 3 | Yes | 14 |
| PF01134 | Glucose inhibited division protein A | 12 | 1 | 13 | 3 | 1 | Yes | 11 |
| PF01209 | ubiE/COQ5 methyltransferase family | 7 | 6 | 13 | 5 | 4 | Yes | 13 |
| PF01546 | Peptidase family M20/M25/M40 | 12 | 1 | 13 | 1 | 1 | No | 0 |
| PF02470 | MlaD protein | 12 | 1 | 13 | 2 | 1 | No | 0 |
| PF02606 | Tetraacyldisaccharide-1-P 4'-kinase | 8 | 5 | 13 | 2 | 1 | No | 0 |
| PF03364 | Polyketide cyclase / dehydrase and lipid transport | 7 | 6 | 13 | 4 | 4 | Yes | 11 |
| PF03848 | Tellurite resistance protein TehB | 6 | 7 | 13 | 2 | 4 | Yes | 13 |
| PF07687 | Peptidase dimerisation domain | 12 | 1 | 13 | 1 | 1 | No | 0 |
| PF10604 | Polyketide cyclase / dehydrase and lipid transport | 10 | 3 | 13 | 4 | 3 | Yes | 8 |
| PF14539 | Domain of unknown function (DUF4442) | 11 | 2 | 13 | 4 | 1 | Yes | 6 |
| PF00743 | Flavin-binding monooxygenase-like | 8 | 4 | 12 | 1 | 1 | Yes | 12 |
| PF01894 | Uncharacterised protein family UPF0047 | 9 | 3 | 12 | 1 | 1 | No | 0 |
| PF03881 | Fructosamine kinase | 6 | 6 | 12 | 1 | 1 | Yes | 12 |
| PF05494 | MlaC protein | 11 | 1 | 12 | 2 | 1 | No | 0 |
| PF08327 | Activator of Hsp90 ATPase homolog 1-like protein | 9 | 3 | 12 | 5 | 2 | No | 0 |
| PF13434 | L-lysine 6-monooxygenase/L-ornithine 5-monooxygenase | 7 | 5 | 12 | 1 | 2 | Yes | 12 |
| PF18693 | TRAM domain | 8 | 4 | 12 | 1 | 1 | Yes | 2 |
| PF01263 | Aldose 1-epimerase | 5 | 6 | 11 | 1 | 1 | No | 0 |
| PF02421 | Ferrous iron transport protein B | 10 | 1 | 11 | 3 | 1 | Yes | 11 |
| PF04909 | Amidohydrolase | 3 | 8 | 11 | 3 | 5 | No | 0 |
| PF01946 | Thi4 family | 4 | 6 | 10 | 3 | 6 | Yes | 10 |
| PF04432 | Coenzyme F420 hydrogenase/dehydrogenase, beta subunit C terminus | 4 | 6 | 10 | 3 | 3 | No | 0 |
| PF12368 | Rhodanase C-terminal | 9 | 1 | 10 | 1 | 1 | No | 0 |
| PF12706 | Beta-lactamase superfamily domain | 5 | 5 | 10 | 3 | 2 | Yes | 7 |
| PF13671 | AAA domain | 8 | 2 | 10 | 4 | 2 | Yes | 6 |
| PF21688 | FAD-dependent protein, C-terminal domain-like | 8 | 2 | 10 | 2 | 2 | No | 0 |
| PF00070 | Pyridine nucleotide-disulphide oxidoreductase | 4 | 5 | 9 | 2 | 5 | Yes | 9 |
| PF01336 | OB-fold nucleic acid binding domain | 8 | 1 | 9 | 3 | 1 | Yes | 1 |
| PF01464 | Transglycosylase SLT domain | 8 | 1 | 9 | 3 | 1 | No | 0 |
| PF03692 | Putative zinc- or iron-chelating domain | 3 | 6 | 9 | 3 | 5 | No | 0 |
| PF04422 | Coenzyme F420 hydrogenase/dehydrogenase, beta subunit N-term | 4 | 5 | 9 | 3 | 2 | No | 0 |
| PF00753 | Metallo-beta-lactamase superfamily | 3 | 5 | 8 | 3 | 2 | Yes | 5 |
| PF01512 | Respiratory-chain NADH dehydrogenase 51 Kd subunit | 2 | 6 | 8 | 1 | 1 | No | 0 |
| PF01844 | HNH endonuclease | 3 | 5 | 8 | 3 | 5 | Yes | 4 |
| PF01863 | YgjP-like, metallopeptidase domain | 5 | 3 | 8 | 2 | 1 | Yes | 1 |
| PF03479 | Plants and Prokaryotes Conserved (PCC) domain | 7 | 1 | 8 | 1 | 1 | No | 0 |
| PF04413 | 3-Deoxy-D-manno-octulosonic-acid transferase (kdotransferase) | 5 | 3 | 8 | 2 | 2 | No | 0 |
| PF06172 | Cupin superfamily (DUF985) | 6 | 2 | 8 | 1 | 1 | No | 0 |
| PF06745 | KaiC | 4 | 4 | 8 | 1 | 2 | Yes | 8 |
| PF10589 | NADH-ubiquinone oxidoreductase-F iron-sulfur binding region | 2 | 6 | 8 | 1 | 1 | No | 0 |
| PF13181 | Tetratricopeptide repeat | 7 | 1 | 8 | 6 | 1 | Yes | 7 |
| PF13481 | AAA domain | 4 | 4 | 8 | 1 | 2 | Yes | 8 |
| PF01695 | IstB-like ATP binding protein | 6 | 1 | 7 | 4 | 1 | Yes | 6 |
| PF01784 | Duf34/NIF3 (NGG1p interacting factor 3) | 6 | 1 | 7 | 2 | 1 | No | 0 |
| PF01833 | IPT/TIG domain | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF01990 | ATP synthase (F/14-kDa) subunit | 5 | 2 | 7 | 3 | 2 | No | 0 |
| PF02586 | SOS response associated peptidase (SRAP) | 4 | 3 | 7 | 2 | 2 | No | 0 |
| PF03746 | LamB/YcsF family | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF03969 | AFG1-like ATPase | 4 | 3 | 7 | 2 | 1 | No | 0 |
| PF04371 | Porphyromonas-type peptidyl-arginine deiminase | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF04472 | Cell division protein SepF | 2 | 5 | 7 | 1 | 2 | No | 0 |
| PF04752 | ChaC-like protein | 2 | 5 | 7 | 1 | 1 | Yes | 1 |
| PF08491 | Squalene epoxidase | 4 | 3 | 7 | 1 | 1 | Yes | 7 |
| PF12804 | MobA-like NTP transferase domain | 4 | 3 | 7 | 4 | 3 | Yes | 3 |
| PF13279 | Thioesterase-like superfamily | 6 | 1 | 7 | 4 | 1 | Yes | 6 |
| PF13431 | Tetratricopeptide repeat | 6 | 1 | 7 | 5 | 1 | Yes | 7 |
| PF13452 | N-terminal half of MaoC dehydratase | 4 | 3 | 7 | 2 | 1 | Yes | 3 |
| PF00025 | ADP-ribosylation factor family | 5 | 1 | 6 | 2 | 1 | Yes | 6 |
| PF00308 | Bacterial DnaA ATPAse domain | 5 | 1 | 6 | 3 | 1 | Yes | 6 |
| PF00903 | Glyoxalase/Bleomycin resistance protein/Dioxygenase superfamily | 5 | 1 | 6 | 3 | 1 | Yes | 6 |
| PF01047 | MarR family | 5 | 1 | 6 | 3 | 1 | Yes | 6 |
| PF01381 | Helix-turn-helix | 5 | 1 | 6 | 4 | 1 | Yes | 4 |
| PF01878 | EVE domain | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF02445 | Quinolinate synthetase A protein | 3 | 3 | 6 | 1 | 1 | No | 0 |
| PF02698 | DUF218 domain | 4 | 2 | 6 | 3 | 1 | No | 0 |
| PF04285 | Protein of unknown function (DUF444) | 5 | 1 | 6 | 1 | 1 | No | 0 |
| PF04468 | PSP1 C-terminal conserved region | 3 | 3 | 6 | 1 | 2 | No | 0 |
| PF06154 | CbeA_antitoxin, type IV, cytoskeleton bundling-enhancing factor A | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF08020 | Protein of unknown function (DUF1706)    | 4 | 2 | 6 | 2 | 2 | Yes | 1 |
| PF12802 | MarR family | 5 | 1 | 6 | 2 | 1 | Yes | 6 |
| PF17146 | PIN domain of ribonuclease | 5 | 1 | 6 | 1 | 1 | No | 0 |
| PF18073 | Rubredoxin metal binding domain | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF00037 | 4Fe-4S binding domain | 2 | 3 | 5 | 2 | 2 | Yes | 3 |
| PF00078 | Reverse transcriptase (RNA-dependent DNA polymerase) | 4 | 1 | 5 | 4 | 1 | No | 0 |
| PF00534 | Glycosyl transferases group 1 | 3 | 2 | 5 | 2 | 2 | Yes | 2 |
| PF03029 | Conserved hypothetical ATP binding protein | 1 | 4 | 5 | 1 | 2 | No | 0 |
| PF04820 | Tryptophan halogenase | 4 | 1 | 5 | 3 | 1 | Yes | 5 |
| PF05673 | Protein of unknown function (DUF815) | 1 | 4 | 5 | 1 | 1 | Yes | 5 |
| PF05991 | YacP-like NYN domain | 4 | 1 | 5 | 1 | 1 | No | 0 |
| PF06032 | Protein of unknown function (DUF917), N-terminal | 3 | 2 | 5 | 2 | 2 | No | 0 |
| PF06039 | Malate:quinone oxidoreductase (Mqo) | 3 | 2 | 5 | 1 | 1 | Yes | 5 |
| PF06228 | Haem utilisation ChuX/HutX | 4 | 1 | 5 | 1 | 1 | Yes | 3 |
| PF06325 | Ribosomal protein L11 methyltransferase (PrmA) | 2 | 3 | 5 | 1 | 3 | Yes | 5 |
| PF07720 | Tetratricopeptide repeat | 4 | 1 | 5 | 1 | 1 | No | 0 |
| PF08450 | SMP-30/Gluconolactonase/LRE-like region | 4 | 1 | 5 | 3 | 1 | No | 0 |
| PF13401 | AAA domain | 2 | 3 | 5 | 2 | 2 | Yes | 5 |
| PF13432 | Tetratricopeptide repeat | 4 | 1 | 5 | 3 | 1 | Yes | 5 |
| PF14464 | Prokaryotic homologs of the JAB domain | 4 | 1 | 5 | 1 | 1 | Yes | 2 |
| PF20906 | Protein of unknown function DUF917, C-terminal | 3 | 2 | 5 | 2 | 2 | No | 0 |
| PF22818 | ApeI dehydratase | 3 | 2 | 5 | 1 | 1 | Yes | 1 |
| PF01663 | Type I phosphodiesterase / nucleotide pyrophosphatase | 2 | 2 | 4 | 1 | 2 | Yes | 3 |
| PF01755 | Glycosyltransferase family 25 (LPS biosynthesis protein) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF01869 | BadF/BadG/BcrA/BcrD ATPase family | 2 | 2 | 4 | 1 | 2 | No | 0 |
| PF01906 | Putative heavy-metal-binding | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF01928 | CYTH domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF01989 | Aconitase X swivel domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF03432 | MobA/VirD2-like, nuclease domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF03747 | ADP-ribosylglycohydrolase | 2 | 2 | 4 | 2 | 2 | No | 0 |
| PF04305 | Protein of unknown function (DUF455) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF04551 | GcpE protein | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF05834 | Lycopene cyclase protein | 2 | 2 | 4 | 1 | 2 | Yes | 4 |
| PF06094 | Gamma-glutamyl cyclotransferase, AIG2-like | 3 | 1 | 4 | 3 | 1 | Yes | 3 |
| PF07977 | FabA-like domain | 3 | 1 | 4 | 3 | 1 | Yes | 1 |
| PF07995 | Glucose / Sorbosone dehydrogenase | 1 | 3 | 4 | 1 | 1 | No | 0 |
| PF10531 | SLBB domain | 1 | 3 | 4 | 1 | 1 | Yes | 3 |
| PF10568 | Outer mitochondrial membrane transport complex protein | 3 | 1 | 4 | 1 | 1 | Yes | 4 |
| PF12681 | Glyoxalase-like domain | 3 | 1 | 4 | 1 | 1 | Yes | 4 |
| PF13191 | AAA ATPase domain | 1 | 3 | 4 | 1 | 2 | Yes | 4 |
| PF13412 | Winged helix-turn-helix DNA-binding | 3 | 1 | 4 | 2 | 1 | Yes | 4 |
| PF13454 | FAD-NAD(P)-binding | 3 | 1 | 4 | 2 | 1 | Yes | 3 |
| PF13541 | Subunit ChlI of Mg-chelatase | 3 | 1 | 4 | 1 | 1 | Yes | 1 |
| PF17171 | Glutathione S-transferase, C-terminal domain | 3 | 1 | 4 | 1 | 1 | Yes | 1 |
| PF17172 | Glutathione S-transferase N-terminal domain  | 3 | 1 | 4 | 1 | 1 | Yes | 4 |
| PF22794 | ZPR1 jelly-roll domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF00884 | Sulfatase | 2 | 1 | 3 | 1 | 1 | Yes | 3 |
| PF00925 | GTP cyclohydrolase II | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01135 | Protein-L-isoaspartate(D-aspartate) O-methyltransferase (PCMT) | 1 | 2 | 3 | 1 | 2 | Yes | 3 |
| PF01343 | Peptidase family S49 | 1 | 2 | 3 | 1 | 1 | Yes | 1 |
| PF01717 | Cobalamin-independent synthase, Catalytic domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01852 | START domain | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF01893 | Uncharacterised protein family UPF0058 | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01960 | ArgJ family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF02190 | ATP-dependent protease La (LON) substrate-binding domain  | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF02401 | LytB protein | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF02541 | Ppx/GppA phosphatase family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF03038 | UL95 family | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF03367 | ZPR1 zinc-finger domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF03807 | NADP oxidoreductase coenzyme F420-dependent | 1 | 2 | 3 | 1 | 2 | Yes | 1 |
| PF04245 | 37-kD nucleoid-associated bacterial protein | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF04670 | Gtr1/RagA G protein conserved region | 2 | 1 | 3 | 2 | 1 | Yes | 2 |
| PF05148 | Hypothetical methyltransferase | 2 | 1 | 3 | 1 | 1 | Yes | 3 |
| PF05362 | Lon protease (S16) C-terminal proteolytic domain | 1 | 2 | 3 | 1 | 2 | Yes | 1 |
| PF05908 | Poly-gamma-glutamate hydrolase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF06080 | Protein of unknown function (DUF938) | 1 | 2 | 3 | 1 | 1 | Yes | 1 |
| PF06627 | Protein of unknown function (DUF1153) | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF07542 | ATP12 chaperone protein | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF07683 | Cobalamin synthesis protein cobW C-terminal domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF08282 | haloacid dehalogenase-like hydrolase | 2 | 1 | 3 | 1 | 1 | Yes | 1 |
| PF08406 | CbbQ/NirQ/NorQ C-terminal  | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF08496 | Peptidase family S49 N-terminal | 1 | 2 | 3 | 1 | 1 | Yes | 1 |
| PF08662 | Eukaryotic translation initiation factor eIF2A | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF09084 | NMT1/THI5 like | 2 | 1 | 3 | 2 | 1 | Yes | 3 |
| PF12917 | 5'-deoxynucleotidase YfbR-like | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF13023 | HD domain | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF13183 | 4Fe-4S dicluster domain | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF13379 | NMT1-like family | 1 | 2 | 3 | 1 | 2 | Yes | 2 |
| PF13460 | NAD(P)H-binding  | 2 | 1 | 3 | 1 | 1 | Yes | 2 |
| PF13483 | Beta-lactamase superfamily domain | 2 | 1 | 3 | 1 | 1 | Yes | 3 |
| PF16347 | N-sulphoglucosamine sulphohydrolase, C-terminal | 2 | 1 | 3 | 1 | 1 | Yes | 3 |
| PF22461 | SLBB domain | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF22817 | ApeP dehydratase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF00201 | UDP-glucoronosyl and UDP-glucosyl transferase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF00436 | Single-strand binding protein family | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF00483 | Nucleotidyl transferase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF00665 | Integrase core domain | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF00702 | haloacid dehalogenase-like hydrolase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF00891 | O-methyltransferase domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF01022 | Bacterial regulatory protein, arsR family | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF01435 | Peptidase family M48 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01564 | Spermine/spermidine synthase domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01774 | UreD urease accessory protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01798 | snoRNA binding domain, fibrillarin | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01866 | Putative diphthamide synthesis protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02550 | Acetyl-CoA hydrolase/transferase N-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02582 | RMND1/Sif2-Sif3/Mrx10, DUF155 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02900 | Catalytic LigB subunit of aromatic ring-opening dioxygenase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02934 | GatB/GatE catalytic domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF03998 | Utp11 protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF04266 | ASCH domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05368 | NmrA-like family | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF05787 | Alkaline phosphatase PhoX | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05812 | Herpesvirus BLRF2 protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05971 | RNA methyltransferase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06224 | DNA glycosylase AlkZ-like | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06722 | Erythromycin biosynthesis protein CIII-like, C-terminal domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF06968 | Biotin and Thiamin Synthesis associated domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07726 | ATPase family associated with various cellular activities (AAA) | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF07978 | NIPSNAP  | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF08388 | Group II intron, maturase-specific domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF12844 | Helix-turn-helix domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13336 | Acetyl-CoA hydrolase/transferase C-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF13391 | HNH endonuclease | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13424 | Tetratricopeptide repeat | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13534 | 4Fe-4S dicluster domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13560 | Helix-turn-helix domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF14697 | 4Fe-4S dicluster domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF22622 | MFE-2 hydratase 2 N-terminal domain    | 1 | 1 | 2 | 1 | 1 | Yes | 2 |

### validation: 81 qualifying families

[Complete CSV](results/split_mixed_label_families/validation_all_GA_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/validation_all_GA_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF13450 | NAD(P)-binding Rossmann-like domain | 8 | 17 | 25 | 4 | 8 | Yes | 24 |
| PF01266 | FAD dependent oxidoreductase | 4 | 17 | 21 | 4 | 7 | Yes | 15 |
| PF07992 | Pyridine nucleotide-disulphide oxidoreductase | 10 | 11 | 21 | 6 | 7 | Yes | 19 |
| PF01494 | FAD binding domain | 5 | 15 | 20 | 3 | 6 | Yes | 20 |
| PF00890 | FAD binding domain | 7 | 11 | 18 | 3 | 7 | Yes | 18 |
| PF00535 | Glycosyl transferase family 2 | 12 | 5 | 17 | 9 | 5 | Yes | 11 |
| PF20143 | ATP-NAD kinase C-terminal domain | 12 | 5 | 17 | 2 | 3 | No | 0 |
| PF00575 | S1 RNA binding domain | 13 | 1 | 14 | 1 | 1 | No | 0 |
| PF01656 | CobQ/CobB/MinD/ParA nucleotide binding domain | 6 | 7 | 13 | 3 | 2 | Yes | 9 |
| PF07728 | AAA domain (dynein-related subfamily) | 9 | 4 | 13 | 1 | 2 | Yes | 11 |
| PF01513 | ATP-NAD kinase N-terminal domain | 8 | 4 | 12 | 2 | 3 | No | 0 |
| PF01926 | 50S ribosome-binding GTPase | 11 | 1 | 12 | 1 | 1 | Yes | 10 |
| PF02622 | Uncharacterized ACR, COG1678 | 11 | 1 | 12 | 1 | 1 | No | 0 |
| PF12831 | FAD dependent oxidoreductase | 5 | 7 | 12 | 3 | 3 | Yes | 12 |
| PF00004 | ATPase family associated with various cellular activities (AAA) | 10 | 1 | 11 | 3 | 1 | Yes | 11 |
| PF01933 | 2-phospho-L-lactate transferase CofD | 8 | 3 | 11 | 1 | 1 | No | 0 |
| PF07726 | ATPase family associated with various cellular activities (AAA) | 9 | 2 | 11 | 1 | 1 | Yes | 11 |
| PF17863 | AAA lid domain | 9 | 2 | 11 | 1 | 1 | No | 0 |
| PF20030 | MoxR domain in the MoxR-vWA-beta-propeller ternary systems | 9 | 2 | 11 | 1 | 1 | Yes | 11 |
| PF07685 | CobB/CobQ-like glutamine amidotransferase domain | 5 | 5 | 10 | 2 | 1 | No | 0 |
| PF13532 | 2OG-Fe(II) oxygenase superfamily | 5 | 5 | 10 | 1 | 2 | No | 0 |
| PF13641 | Glycosyltransferase like family 2 | 8 | 2 | 10 | 6 | 2 | Yes | 10 |
| PF00005 | ABC transporter | 5 | 4 | 9 | 2 | 1 | Yes | 2 |
| PF00472 | RF-1 domain | 8 | 1 | 9 | 1 | 1 | No | 0 |
| PF03577 | Peptidase family C69 | 5 | 4 | 9 | 1 | 1 | No | 0 |
| PF13738 | Pyridine nucleotide-disulphide oxidoreductase | 5 | 4 | 9 | 4 | 2 | Yes | 9 |
| PF13806 | Rieske-like [2Fe-2S] domain | 6 | 3 | 9 | 1 | 1 | Yes | 3 |
| PF01134 | Glucose inhibited division protein A | 3 | 5 | 8 | 2 | 3 | Yes | 8 |
| PF02492 | CobW/HypB/UreG, nucleotide-binding domain | 4 | 4 | 8 | 2 | 2 | No | 0 |
| PF13489 | Methyltransferase domain | 1 | 6 | 7 | 1 | 2 | Yes | 7 |
| PF01593 | Flavin containing amine oxidoreductase | 3 | 3 | 6 | 1 | 3 | Yes | 6 |
| PF03008 | Archaea bacterial proteins of unknown function | 2 | 4 | 6 | 2 | 1 | No | 0 |
| PF05257 | CHAP domain | 4 | 2 | 6 | 3 | 1 | No | 0 |
| PF13181 | Tetratricopeptide repeat | 5 | 1 | 6 | 4 | 1 | Yes | 6 |
| PF13500 | AAA domain | 4 | 2 | 6 | 2 | 1 | Yes | 6 |
| PF01637 | ATPase domain predominantly from Archaea | 1 | 4 | 5 | 1 | 1 | Yes | 5 |
| PF01813 | ATP synthase subunit D  | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF04363 | Protein of unknown function (DUF496) | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF05834 | Lycopene cyclase protein | 3 | 2 | 5 | 1 | 2 | Yes | 5 |
| PF06838 | Methionine gamma-lyase  | 4 | 1 | 5 | 1 | 1 | Yes | 4 |
| PF08327 | Activator of Hsp90 ATPase homolog 1-like protein | 4 | 1 | 5 | 2 | 1 | Yes | 2 |
| PF13173 | AAA domain | 1 | 4 | 5 | 1 | 1 | Yes | 5 |
| PF13434 | L-lysine 6-monooxygenase/L-ornithine 5-monooxygenase | 3 | 2 | 5 | 2 | 1 | Yes | 5 |
| PF14559 | Tetratricopeptide repeat | 4 | 1 | 5 | 4 | 1 | Yes | 4 |
| PF01053 | Cys/Met metabolism PLP-dependent enzyme | 3 | 1 | 4 | 1 | 1 | Yes | 4 |
| PF01476 | LysM domain | 2 | 2 | 4 | 2 | 1 | No | 0 |
| PF01844 | HNH endonuclease | 3 | 1 | 4 | 3 | 1 | Yes | 2 |
| PF04536 | TPM domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF04991 | LicD family | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF14691 | Dihydroprymidine dehydrogenase domain II, 4Fe-4S cluster | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF22578 | Geranylgeranyl reductase catalytic domain | 2 | 2 | 4 | 1 | 1 | Yes | 3 |
| PF00753 | Metallo-beta-lactamase superfamily | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF01139 | tRNA-splicing ligase RtcB | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01519 | Protein of unknown function DUF16 | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF02374 | Anion-transporting ATPase | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF03437 | BtpA family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF04055 | Radical SAM superfamily | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF08014 | Microtubule-associated tyrosine carboxypeptidase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF09350 | DnaJ homologue, subfamily C, member 28, conserved domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF13401 | AAA domain | 2 | 1 | 3 | 2 | 1 | Yes | 3 |
| PF13454 | FAD-NAD(P)-binding | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF13614 | AAA domain | 1 | 2 | 3 | 1 | 1 | Yes | 3 |
| PF00340 | Interleukin-1 / 18 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01435 | Peptidase family M48 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02709 | N-terminal domain of galactosyltransferase | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF03364 | Polyketide cyclase / dehydrase and lipid transport | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF04280 | Tim44-like domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05240 | APOBEC-like C-terminal domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF07683 | Cobalamin synthesis protein cobW C-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF08210 | APOBEC-like N-terminal domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF12895 | Anaphase-promoting complex, cyclosome, subunit 3 | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF13428 | Tetratricopeptide repeat | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF14279 | HNH endonuclease | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18750 | Secreted Novel AID/APOBEC-like Deaminase 4 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18769 | APOBEC1 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18771 | APOBEC3 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18772 | APOBEC2 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18774 | APOBEC4-like -AID/APOBEC-deaminase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18775 | APOBEC4 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18778 | Novel AID APOBEC clade 1 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF18782 | Novel AID APOBEC clade 2 | 1 | 1 | 2 | 1 | 1 | Yes | 2 |

### reserved_test: 100 qualifying families

[Complete CSV](results/split_mixed_label_families/reserved_test_all_GA_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/reserved_test_all_GA_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF01494 | FAD binding domain | 11 | 16 | 27 | 2 | 3 | Yes | 18 |
| PF13450 | NAD(P)-binding Rossmann-like domain | 7 | 18 | 25 | 3 | 6 | Yes | 25 |
| PF01527 | Transposase | 21 | 1 | 22 | 1 | 1 | Yes | 8 |
| PF13460 | NAD(P)H-binding  | 17 | 5 | 22 | 1 | 1 | Yes | 17 |
| PF13304 | AAA domain, putative AbiEii toxin, Type IV TA system | 6 | 13 | 19 | 4 | 3 | Yes | 12 |
| PF00890 | FAD binding domain | 7 | 11 | 18 | 4 | 7 | Yes | 17 |
| PF01370 | NAD dependent epimerase/dehydratase family | 15 | 1 | 16 | 1 | 1 | Yes | 16 |
| PF01715 | IPP transferase | 14 | 2 | 16 | 1 | 1 | Yes | 14 |
| PF13476 | AAA domain | 3 | 12 | 15 | 2 | 3 | Yes | 10 |
| PF01745 | Isopentenyl transferase | 13 | 1 | 14 | 1 | 1 | Yes | 14 |
| PF05368 | NmrA-like family | 12 | 2 | 14 | 1 | 1 | Yes | 14 |
| PF00534 | Glycosyl transferases group 1 | 9 | 4 | 13 | 3 | 3 | Yes | 13 |
| PF12804 | MobA-like NTP transferase domain | 6 | 7 | 13 | 3 | 3 | Yes | 6 |
| PF13692 | Glycosyl transferases group 1 | 9 | 4 | 13 | 3 | 3 | Yes | 13 |
| PF00535 | Glycosyl transferase family 2 | 10 | 2 | 12 | 7 | 2 | Yes | 10 |
| PF01266 | FAD dependent oxidoreductase | 3 | 9 | 12 | 2 | 4 | Yes | 12 |
| PF00753 | Metallo-beta-lactamase superfamily | 8 | 3 | 11 | 4 | 3 | Yes | 7 |
| PF03618 | Kinase/pyrophosphorylase | 8 | 3 | 11 | 1 | 1 | No | 0 |
| PF01977 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase Rift-related domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF13439 | Glycosyltransferase Family 4 | 7 | 3 | 10 | 2 | 2 | Yes | 10 |
| PF13579 | Glycosyl transferase 4-like domain | 7 | 3 | 10 | 2 | 2 | Yes | 10 |
| PF13641 | Glycosyltransferase like family 2 | 8 | 2 | 10 | 5 | 2 | Yes | 10 |
| PF20695 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase N-terminal domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF20696 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase C-terminal domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF02410 | Ribosomal silencing factor during starvation  | 6 | 3 | 9 | 1 | 2 | No | 0 |
| PF03054 | tRNA methyl transferase HUP domain | 7 | 2 | 9 | 1 | 1 | Yes | 8 |
| PF20258 | Aminomethyltransferase beta-barrel domain | 7 | 2 | 9 | 1 | 1 | No | 0 |
| PF20259 | tRNA methyl transferase PRC-barrel domain | 7 | 2 | 9 | 1 | 1 | No | 0 |
| PF01695 | IstB-like ATP binding protein | 3 | 5 | 8 | 1 | 3 | Yes | 6 |
| PF01812 | 5-formyltetrahydrofolate cyclo-ligase family | 6 | 2 | 8 | 2 | 1 | No | 0 |
| PF12706 | Beta-lactamase superfamily domain | 6 | 2 | 8 | 3 | 2 | Yes | 7 |
| PF01435 | Peptidase family M48 | 4 | 3 | 7 | 1 | 2 | No | 0 |
| PF02540 | NAD synthase | 5 | 2 | 7 | 1 | 1 | Yes | 7 |
| PF03942 | DTW domain | 3 | 4 | 7 | 1 | 1 | No | 0 |
| PF00149 | Calcineurin-like phosphoesterase | 1 | 5 | 6 | 1 | 3 | Yes | 3 |
| PF00308 | Bacterial DnaA ATPAse domain | 5 | 1 | 6 | 2 | 1 | Yes | 6 |
| PF01774 | UreD urease accessory protein | 1 | 5 | 6 | 1 | 2 | No | 0 |
| PF08241 | Methyltransferase domain | 3 | 3 | 6 | 3 | 3 | Yes | 5 |
| PF10111 | Glycosyltransferase like family 2 | 5 | 1 | 6 | 2 | 1 | Yes | 6 |
| PF12831 | FAD dependent oxidoreductase | 3 | 3 | 6 | 2 | 3 | Yes | 6 |
| PF13175 | AAA ATPase domain | 2 | 4 | 6 | 2 | 3 | Yes | 6 |
| PF20706 | Family 4 Glycosyltransferase in conflict systems | 3 | 3 | 6 | 3 | 2 | Yes | 6 |
| PF00004 | ATPase family associated with various cellular activities (AAA) | 1 | 4 | 5 | 1 | 3 | Yes | 5 |
| PF02626 | Carboxyltransferase domain, subdomain A and B  | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF03435 | Saccharopine dehydrogenase NADP binding domain | 4 | 1 | 5 | 2 | 1 | Yes | 2 |
| PF08242 | Methyltransferase domain | 3 | 2 | 5 | 3 | 2 | Yes | 5 |
| PF10604 | Polyketide cyclase / dehydrase and lipid transport | 3 | 2 | 5 | 1 | 2 | Yes | 4 |
| PF13181 | Tetratricopeptide repeat | 4 | 1 | 5 | 2 | 1 | Yes | 5 |
| PF13432 | Tetratricopeptide repeat | 4 | 1 | 5 | 3 | 1 | Yes | 5 |
| PF13477 | Glycosyl transferase 4-like | 4 | 1 | 5 | 2 | 1 | Yes | 5 |
| PF13649 | Methyltransferase domain | 3 | 2 | 5 | 3 | 2 | Yes | 5 |
| PF13847 | Methyltransferase domain | 4 | 1 | 5 | 4 | 1 | Yes | 5 |
| PF01676 | Metalloenzyme superfamily | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF02811 | PHP domain | 3 | 1 | 4 | 2 | 1 | No | 0 |
| PF03061 | Thioesterase superfamily | 1 | 3 | 4 | 1 | 3 | Yes | 2 |
| PF03358 | NADPH-dependent FMN reductase | 2 | 2 | 4 | 1 | 1 | Yes | 2 |
| PF03364 | Polyketide cyclase / dehydrase and lipid transport | 3 | 1 | 4 | 1 | 1 | Yes | 4 |
| PF07081 | Protein of unknown function (DUF1349) | 3 | 1 | 4 | 1 | 1 | Yes | 2 |
| PF07451 | Stage V sporulation protein AD (SpoVAD) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF09084 | NMT1/THI5 like | 3 | 1 | 4 | 3 | 1 | Yes | 2 |
| PF10143 | 2,3-bisphosphoglycerate-independent phosphoglycerate mutase | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF10263 | SprT-like family | 2 | 2 | 4 | 2 | 1 | No | 0 |
| PF00037 | 4Fe-4S binding domain | 1 | 2 | 3 | 1 | 2 | No | 0 |
| PF00355 | Rieske [2Fe-2S] domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01128 | 2-C-methyl-D-erythritol 4-phosphate cytidylyltransferase | 1 | 2 | 3 | 1 | 2 | Yes | 3 |
| PF01142 | tRNA pseudouridine synthase D (TruD) | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01343 | Peptidase family S49 | 2 | 1 | 3 | 1 | 1 | Yes | 2 |
| PF01436 | NHL repeat | 1 | 2 | 3 | 1 | 2 | No | 0 |
| PF02636 | Putative S-adenosyl-L-methionine-dependent methyltransferase | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF07676 | WD40-like Beta Propeller Repeat | 2 | 1 | 3 | 2 | 1 | Yes | 1 |
| PF07719 | Tetratricopeptide repeat | 2 | 1 | 3 | 2 | 1 | Yes | 3 |
| PF10275 | Peptidase C65 Otubain | 1 | 2 | 3 | 1 | 1 | Yes | 2 |
| PF13489 | Methyltransferase domain | 2 | 1 | 3 | 2 | 1 | Yes | 3 |
| PF13524 | Glycosyl transferases group 1 | 2 | 1 | 3 | 2 | 1 | Yes | 3 |
| PF16653 | Saccharopine dehydrogenase C-terminal domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF17784 | Sulfotransferase domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF00069 | Protein kinase domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF00515 | Tetratricopeptide repeat | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF00574 | Clp protease | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF01656 | CobQ/CobB/MinD/ParA nucleotide binding domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF01734 | Patatin-like phospholipase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01809 | Putative membrane protein insertion efficiency factor | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01968 | Hydantoinase/oxoprolinase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02374 | Anion-transporting ATPase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF02470 | MlaD protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02581 | Thiamine monophosphate synthase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF04909 | Amidohydrolase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05257 | CHAP domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05378 | Hydantoinase/oxoprolinase N-terminal region | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06089 | L-asparaginase II | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07714 | Protein tyrosine and serine/threonine kinase | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF12850 | Calcineurin-like phosphoesterase superfamily domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13174 | Tetratricopeptide repeat | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13263 | PHP-associated | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF13279 | Thioesterase-like superfamily | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF13414 | TPR repeat | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13519 | von Willebrand factor type A domain | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF13614 | AAA domain | 1 | 1 | 2 | 1 | 1 | Yes | 2 |
| PF13671 | AAA domain | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF17283 | SprT-like zinc ribbon domain | 1 | 1 | 2 | 1 | 1 | No | 0 |


## Sensitivity: remove competing overlapping hits

### train: 150 qualifying families

[Complete CSV](results/split_mixed_label_families/train_no_competing_overlap_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/train_no_competing_overlap_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF00455 | DeoR C terminal sensor domain | 31 | 26 | 57 | 1 | 1 | No | 0 |
| PF04055 | Radical SAM superfamily | 31 | 21 | 52 | 4 | 11 | No | 0 |
| PF00575 | S1 RNA binding domain | 48 | 2 | 50 | 6 | 1 | Yes | 2 |
| PF00919 | Uncharacterized protein family UPF0004 | 25 | 11 | 36 | 1 | 1 | No | 0 |
| PF07282 | Putative transposase DNA-binding domain | 22 | 6 | 28 | 10 | 5 | No | 0 |
| PF01966 | HD domain | 22 | 2 | 24 | 2 | 1 | No | 0 |
| PF01588 | Putative tRNA binding domain | 20 | 2 | 22 | 2 | 1 | No | 0 |
| PF01938 | TRAM domain | 15 | 7 | 22 | 2 | 1 | No | 0 |
| PF01368 | DHH family | 16 | 2 | 18 | 1 | 1 | No | 0 |
| PF01385 | Probable transposase | 17 | 1 | 18 | 5 | 1 | No | 0 |
| PF13167 | GTP-binding GTPase N-terminal | 14 | 4 | 18 | 1 | 1 | No | 0 |
| PF16360 | GTP-binding GTPase Middle Region | 14 | 4 | 18 | 1 | 1 | No | 0 |
| PF02272 | DHHA1 domain | 14 | 2 | 16 | 1 | 1 | No | 0 |
| PF03966 | Trm112p-like protein | 15 | 1 | 16 | 2 | 1 | No | 0 |
| PF04390 | Lipopolysaccharide-assembly | 15 | 1 | 16 | 2 | 1 | No | 0 |
| PF08220 | DeoR-like helix-turn-helix domain | 7 | 9 | 16 | 1 | 1 | No | 0 |
| PF00535 | Glycosyl transferase family 2 | 13 | 2 | 15 | 8 | 2 | No | 0 |
| PF00696 | Amino acid kinase family | 11 | 4 | 15 | 2 | 3 | No | 0 |
| PF01266 | FAD dependent oxidoreductase | 9 | 6 | 15 | 2 | 2 | No | 0 |
| PF16901 | C-terminal domain of alpha-glycerophosphate oxidase | 10 | 5 | 15 | 1 | 1 | No | 0 |
| PF01636 | Phosphotransferase enzyme family | 8 | 6 | 14 | 3 | 4 | No | 0 |
| PF02666 | Phosphatidylserine decarboxylase | 9 | 5 | 14 | 2 | 3 | No | 0 |
| PF04012 | PspA/IM30 family | 7 | 7 | 14 | 1 | 1 | No | 0 |
| PF13353 | 4Fe-4S single cluster domain | 10 | 4 | 14 | 2 | 2 | No | 0 |
| PF01546 | Peptidase family M20/M25/M40 | 12 | 1 | 13 | 1 | 1 | No | 0 |
| PF02470 | MlaD protein | 12 | 1 | 13 | 2 | 1 | No | 0 |
| PF02606 | Tetraacyldisaccharide-1-P 4'-kinase | 8 | 5 | 13 | 2 | 1 | No | 0 |
| PF07687 | Peptidase dimerisation domain | 12 | 1 | 13 | 1 | 1 | No | 0 |
| PF07728 | AAA domain (dynein-related subfamily) | 4 | 9 | 13 | 3 | 5 | No | 0 |
| PF01494 | FAD binding domain | 3 | 9 | 12 | 1 | 5 | Yes | 1 |
| PF01894 | Uncharacterised protein family UPF0047 | 9 | 3 | 12 | 1 | 1 | No | 0 |
| PF05494 | MlaC protein | 11 | 1 | 12 | 2 | 1 | No | 0 |
| PF08327 | Activator of Hsp90 ATPase homolog 1-like protein | 9 | 3 | 12 | 5 | 2 | No | 0 |
| PF01263 | Aldose 1-epimerase | 5 | 6 | 11 | 1 | 1 | No | 0 |
| PF04909 | Amidohydrolase | 3 | 8 | 11 | 3 | 5 | No | 0 |
| PF04432 | Coenzyme F420 hydrogenase/dehydrogenase, beta subunit C terminus | 4 | 6 | 10 | 3 | 3 | No | 0 |
| PF12368 | Rhodanase C-terminal | 9 | 1 | 10 | 1 | 1 | No | 0 |
| PF18693 | TRAM domain | 6 | 4 | 10 | 1 | 1 | No | 0 |
| PF21688 | FAD-dependent protein, C-terminal domain-like | 8 | 2 | 10 | 2 | 2 | No | 0 |
| PF01464 | Transglycosylase SLT domain | 8 | 1 | 9 | 3 | 1 | No | 0 |
| PF03692 | Putative zinc- or iron-chelating domain | 3 | 6 | 9 | 3 | 5 | No | 0 |
| PF04422 | Coenzyme F420 hydrogenase/dehydrogenase, beta subunit N-term | 4 | 5 | 9 | 3 | 2 | No | 0 |
| PF01336 | OB-fold nucleic acid binding domain | 7 | 1 | 8 | 3 | 1 | No | 0 |
| PF01512 | Respiratory-chain NADH dehydrogenase 51 Kd subunit | 2 | 6 | 8 | 1 | 1 | No | 0 |
| PF01926 | 50S ribosome-binding GTPase | 6 | 2 | 8 | 1 | 1 | No | 0 |
| PF03479 | Plants and Prokaryotes Conserved (PCC) domain | 7 | 1 | 8 | 1 | 1 | No | 0 |
| PF04413 | 3-Deoxy-D-manno-octulosonic-acid transferase (kdotransferase) | 5 | 3 | 8 | 2 | 2 | No | 0 |
| PF06172 | Cupin superfamily (DUF985) | 6 | 2 | 8 | 1 | 1 | No | 0 |
| PF10589 | NADH-ubiquinone oxidoreductase-F iron-sulfur binding region | 2 | 6 | 8 | 1 | 1 | No | 0 |
| PF01784 | Duf34/NIF3 (NGG1p interacting factor 3) | 6 | 1 | 7 | 2 | 1 | No | 0 |
| PF01833 | IPT/TIG domain | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF01863 | YgjP-like, metallopeptidase domain | 4 | 3 | 7 | 2 | 1 | No | 0 |
| PF01990 | ATP synthase (F/14-kDa) subunit | 5 | 2 | 7 | 3 | 2 | No | 0 |
| PF02586 | SOS response associated peptidase (SRAP) | 4 | 3 | 7 | 2 | 2 | No | 0 |
| PF03746 | LamB/YcsF family | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF03969 | AFG1-like ATPase | 4 | 3 | 7 | 2 | 1 | No | 0 |
| PF04371 | Porphyromonas-type peptidyl-arginine deiminase | 6 | 1 | 7 | 1 | 1 | No | 0 |
| PF04472 | Cell division protein SepF | 2 | 5 | 7 | 1 | 2 | No | 0 |
| PF14539 | Domain of unknown function (DUF4442) | 5 | 2 | 7 | 2 | 1 | No | 0 |
| PF01878 | EVE domain | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF02445 | Quinolinate synthetase A protein | 3 | 3 | 6 | 1 | 1 | No | 0 |
| PF02698 | DUF218 domain | 4 | 2 | 6 | 3 | 1 | No | 0 |
| PF03061 | Thioesterase superfamily | 5 | 1 | 6 | 2 | 1 | No | 0 |
| PF04285 | Protein of unknown function (DUF444) | 5 | 1 | 6 | 1 | 1 | No | 0 |
| PF04468 | PSP1 C-terminal conserved region | 3 | 3 | 6 | 1 | 2 | No | 0 |
| PF04752 | ChaC-like protein | 1 | 5 | 6 | 1 | 1 | No | 0 |
| PF06154 | CbeA_antitoxin, type IV, cytoskeleton bundling-enhancing factor A | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF07992 | Pyridine nucleotide-disulphide oxidoreductase | 2 | 4 | 6 | 2 | 3 | No | 0 |
| PF17146 | PIN domain of ribonuclease | 5 | 1 | 6 | 1 | 1 | No | 0 |
| PF18073 | Rubredoxin metal binding domain | 4 | 2 | 6 | 1 | 1 | No | 0 |
| PF00078 | Reverse transcriptase (RNA-dependent DNA polymerase) | 4 | 1 | 5 | 4 | 1 | No | 0 |
| PF03029 | Conserved hypothetical ATP binding protein | 1 | 4 | 5 | 1 | 2 | No | 0 |
| PF05991 | YacP-like NYN domain | 4 | 1 | 5 | 1 | 1 | No | 0 |
| PF06032 | Protein of unknown function (DUF917), N-terminal | 3 | 2 | 5 | 2 | 2 | No | 0 |
| PF07720 | Tetratricopeptide repeat | 4 | 1 | 5 | 1 | 1 | No | 0 |
| PF08020 | Protein of unknown function (DUF1706)    | 3 | 2 | 5 | 1 | 2 | No | 0 |
| PF08450 | SMP-30/Gluconolactonase/LRE-like region | 4 | 1 | 5 | 3 | 1 | No | 0 |
| PF10604 | Polyketide cyclase / dehydrase and lipid transport | 4 | 1 | 5 | 1 | 1 | No | 0 |
| PF20906 | Protein of unknown function DUF917, C-terminal | 3 | 2 | 5 | 2 | 2 | No | 0 |
| PF01755 | Glycosyltransferase family 25 (LPS biosynthesis protein) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF01844 | HNH endonuclease | 1 | 3 | 4 | 1 | 3 | No | 0 |
| PF01869 | BadF/BadG/BcrA/BcrD ATPase family | 2 | 2 | 4 | 1 | 2 | No | 0 |
| PF01906 | Putative heavy-metal-binding | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF01928 | CYTH domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF01989 | Aconitase X swivel domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF02374 | Anion-transporting ATPase | 1 | 3 | 4 | 1 | 1 | No | 0 |
| PF03432 | MobA/VirD2-like, nuclease domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF03747 | ADP-ribosylglycohydrolase | 2 | 2 | 4 | 2 | 2 | No | 0 |
| PF04305 | Protein of unknown function (DUF455) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF04551 | GcpE protein | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF07995 | Glucose / Sorbosone dehydrogenase | 1 | 3 | 4 | 1 | 1 | No | 0 |
| PF12804 | MobA-like NTP transferase domain | 2 | 2 | 4 | 2 | 2 | No | 0 |
| PF13452 | N-terminal half of MaoC dehydratase | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF13671 | AAA domain | 2 | 2 | 4 | 1 | 2 | No | 0 |
| PF22794 | ZPR1 jelly-roll domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF22818 | ApeI dehydratase | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF00534 | Glycosyl transferases group 1 | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF00753 | Metallo-beta-lactamase superfamily | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF00925 | GTP cyclohydrolase II | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01656 | CobQ/CobB/MinD/ParA nucleotide binding domain | 1 | 2 | 3 | 1 | 2 | No | 0 |
| PF01717 | Cobalamin-independent synthase, Catalytic domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01893 | Uncharacterised protein family UPF0058 | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01960 | ArgJ family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF02190 | ATP-dependent protease La (LON) substrate-binding domain  | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF02401 | LytB protein | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF02541 | Ppx/GppA phosphatase family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF03038 | UL95 family | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF03367 | ZPR1 zinc-finger domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF04245 | 37-kD nucleoid-associated bacterial protein | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF05908 | Poly-gamma-glutamate hydrolase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF06627 | Protein of unknown function (DUF1153) | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF07542 | ATP12 chaperone protein | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF07683 | Cobalamin synthesis protein cobW C-terminal domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF07977 | FabA-like domain | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF08406 | CbbQ/NirQ/NorQ C-terminal  | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF08662 | Eukaryotic translation initiation factor eIF2A | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF12706 | Beta-lactamase superfamily domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF14464 | Prokaryotic homologs of the JAB domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF17171 | Glutathione S-transferase, C-terminal domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF22817 | ApeP dehydratase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF00037 | 4Fe-4S binding domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF00436 | Single-strand binding protein family | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01343 | Peptidase family S49 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01435 | Peptidase family M48 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01564 | Spermine/spermidine synthase domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01774 | UreD urease accessory protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01798 | snoRNA binding domain, fibrillarin | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01866 | Putative diphthamide synthesis protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02550 | Acetyl-CoA hydrolase/transferase N-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02582 | RMND1/Sif2-Sif3/Mrx10, DUF155 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02900 | Catalytic LigB subunit of aromatic ring-opening dioxygenase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02934 | GatB/GatE catalytic domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF03807 | NADP oxidoreductase coenzyme F420-dependent | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF03998 | Utp11 protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF04266 | ASCH domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05362 | Lon protease (S16) C-terminal proteolytic domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05787 | Alkaline phosphatase PhoX | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05812 | Herpesvirus BLRF2 protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05971 | RNA methyltransferase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06080 | Protein of unknown function (DUF938) | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06094 | Gamma-glutamyl cyclotransferase, AIG2-like | 1 | 1 | 2 | 1 | 1 | Yes | 1 |
| PF06224 | DNA glycosylase AlkZ-like | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06228 | Haem utilisation ChuX/HutX | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06968 | Biotin and Thiamin Synthesis associated domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07978 | NIPSNAP  | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF08282 | haloacid dehalogenase-like hydrolase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF08388 | Group II intron, maturase-specific domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF08496 | Peptidase family S49 N-terminal | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF13336 | Acetyl-CoA hydrolase/transferase C-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF13450 | NAD(P)-binding Rossmann-like domain | 1 | 1 | 2 | 1 | 1 | No | 0 |

### validation: 38 qualifying families

[Complete CSV](results/split_mixed_label_families/validation_no_competing_overlap_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/validation_no_competing_overlap_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF20143 | ATP-NAD kinase C-terminal domain | 12 | 5 | 17 | 2 | 3 | No | 0 |
| PF00575 | S1 RNA binding domain | 13 | 1 | 14 | 1 | 1 | No | 0 |
| PF01513 | ATP-NAD kinase N-terminal domain | 8 | 4 | 12 | 2 | 3 | No | 0 |
| PF02622 | Uncharacterized ACR, COG1678 | 11 | 1 | 12 | 1 | 1 | No | 0 |
| PF01933 | 2-phospho-L-lactate transferase CofD | 8 | 3 | 11 | 1 | 1 | No | 0 |
| PF17863 | AAA lid domain | 9 | 2 | 11 | 1 | 1 | No | 0 |
| PF07685 | CobB/CobQ-like glutamine amidotransferase domain | 5 | 5 | 10 | 2 | 1 | No | 0 |
| PF13532 | 2OG-Fe(II) oxygenase superfamily | 5 | 5 | 10 | 1 | 2 | No | 0 |
| PF00472 | RF-1 domain | 8 | 1 | 9 | 1 | 1 | No | 0 |
| PF03577 | Peptidase family C69 | 5 | 4 | 9 | 1 | 1 | No | 0 |
| PF02492 | CobW/HypB/UreG, nucleotide-binding domain | 4 | 4 | 8 | 2 | 2 | No | 0 |
| PF00005 | ABC transporter | 3 | 4 | 7 | 2 | 1 | No | 0 |
| PF00535 | Glycosyl transferase family 2 | 3 | 3 | 6 | 3 | 3 | No | 0 |
| PF01266 | FAD dependent oxidoreductase | 2 | 4 | 6 | 2 | 2 | No | 0 |
| PF03008 | Archaea bacterial proteins of unknown function | 2 | 4 | 6 | 2 | 1 | No | 0 |
| PF05257 | CHAP domain | 4 | 2 | 6 | 3 | 1 | No | 0 |
| PF13806 | Rieske-like [2Fe-2S] domain | 3 | 3 | 6 | 1 | 1 | No | 0 |
| PF01813 | ATP synthase subunit D  | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF04363 | Protein of unknown function (DUF496) | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF01476 | LysM domain | 2 | 2 | 4 | 2 | 1 | No | 0 |
| PF01656 | CobQ/CobB/MinD/ParA nucleotide binding domain | 1 | 3 | 4 | 1 | 1 | No | 0 |
| PF04536 | TPM domain | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF04991 | LicD family | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF14691 | Dihydroprymidine dehydrogenase domain II, 4Fe-4S cluster | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF00753 | Metallo-beta-lactamase superfamily | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF01139 | tRNA-splicing ligase RtcB | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01519 | Protein of unknown function DUF16 | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF03437 | BtpA family | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF04055 | Radical SAM superfamily | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF07992 | Pyridine nucleotide-disulphide oxidoreductase | 1 | 2 | 3 | 1 | 2 | Yes | 1 |
| PF08014 | Microtubule-associated tyrosine carboxypeptidase | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF08327 | Activator of Hsp90 ATPase homolog 1-like protein | 2 | 1 | 3 | 2 | 1 | No | 0 |
| PF09350 | DnaJ homologue, subfamily C, member 28, conserved domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF00340 | Interleukin-1 / 18 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01435 | Peptidase family M48 | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01926 | 50S ribosome-binding GTPase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF04280 | Tim44-like domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07683 | Cobalamin synthesis protein cobW C-terminal domain | 1 | 1 | 2 | 1 | 1 | No | 0 |

### reserved_test: 48 qualifying families

[Complete CSV](results/split_mixed_label_families/reserved_test_no_competing_overlap_mixed_families.csv) · [Complete TSV](results/split_mixed_label_families/reserved_test_no_competing_overlap_mixed_families.tsv)

| Pfam | Biological family/domain name | Soluble | Insoluble | Total | Groups: soluble | Groups: insoluble | Competing overlap* | Affected proteins* |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |
| PF01527 | Transposase | 13 | 1 | 14 | 1 | 1 | No | 0 |
| PF03618 | Kinase/pyrophosphorylase | 8 | 3 | 11 | 1 | 1 | No | 0 |
| PF01494 | FAD binding domain | 5 | 5 | 10 | 1 | 1 | Yes | 1 |
| PF01977 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase Rift-related domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF20695 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase N-terminal domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF20696 | 3-octaprenyl-4-hydroxybenzoate carboxy-lyase C-terminal domain | 2 | 8 | 10 | 1 | 1 | No | 0 |
| PF02410 | Ribosomal silencing factor during starvation  | 6 | 3 | 9 | 1 | 2 | No | 0 |
| PF20258 | Aminomethyltransferase beta-barrel domain | 7 | 2 | 9 | 1 | 1 | No | 0 |
| PF20259 | tRNA methyl transferase PRC-barrel domain | 7 | 2 | 9 | 1 | 1 | No | 0 |
| PF01812 | 5-formyltetrahydrofolate cyclo-ligase family | 6 | 2 | 8 | 2 | 1 | No | 0 |
| PF01435 | Peptidase family M48 | 4 | 3 | 7 | 1 | 2 | No | 0 |
| PF03942 | DTW domain | 3 | 4 | 7 | 1 | 1 | No | 0 |
| PF12804 | MobA-like NTP transferase domain | 4 | 3 | 7 | 2 | 1 | No | 0 |
| PF13304 | AAA domain, putative AbiEii toxin, Type IV TA system | 4 | 3 | 7 | 2 | 1 | No | 0 |
| PF01774 | UreD urease accessory protein | 1 | 5 | 6 | 1 | 2 | No | 0 |
| PF02626 | Carboxyltransferase domain, subdomain A and B  | 1 | 4 | 5 | 1 | 1 | No | 0 |
| PF13460 | NAD(P)H-binding  | 2 | 3 | 5 | 1 | 1 | No | 0 |
| PF13476 | AAA domain | 2 | 3 | 5 | 1 | 1 | No | 0 |
| PF00753 | Metallo-beta-lactamase superfamily | 2 | 2 | 4 | 2 | 2 | No | 0 |
| PF01676 | Metalloenzyme superfamily | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF02811 | PHP domain | 3 | 1 | 4 | 2 | 1 | No | 0 |
| PF07451 | Stage V sporulation protein AD (SpoVAD) | 2 | 2 | 4 | 1 | 1 | No | 0 |
| PF10143 | 2,3-bisphosphoglycerate-independent phosphoglycerate mutase | 3 | 1 | 4 | 1 | 1 | No | 0 |
| PF10263 | SprT-like family | 2 | 2 | 4 | 2 | 1 | No | 0 |
| PF00037 | 4Fe-4S binding domain | 1 | 2 | 3 | 1 | 2 | No | 0 |
| PF00355 | Rieske [2Fe-2S] domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01142 | tRNA pseudouridine synthase D (TruD) | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01436 | NHL repeat | 1 | 2 | 3 | 1 | 2 | No | 0 |
| PF02636 | Putative S-adenosyl-L-methionine-dependent methyltransferase | 1 | 2 | 3 | 1 | 1 | No | 0 |
| PF03435 | Saccharopine dehydrogenase NADP binding domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF16653 | Saccharopine dehydrogenase C-terminal domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF17784 | Sulfotransferase domain | 2 | 1 | 3 | 1 | 1 | No | 0 |
| PF01715 | IPP transferase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01734 | Patatin-like phospholipase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01809 | Putative membrane protein insertion efficiency factor | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF01968 | Hydantoinase/oxoprolinase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02470 | MlaD protein | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF02581 | Thiamine monophosphate synthase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF03061 | Thioesterase superfamily | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF04909 | Amidohydrolase | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05257 | CHAP domain | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF05378 | Hydantoinase/oxoprolinase N-terminal region | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF06089 | L-asparaginase II | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07081 | Protein of unknown function (DUF1349) | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF07676 | WD40-like Beta Propeller Repeat | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF09084 | NMT1/THI5 like | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF13263 | PHP-associated | 1 | 1 | 2 | 1 | 1 | No | 0 |
| PF17283 | SprT-like zinc ribbon domain | 1 | 1 | 2 | 1 | 1 | No | 0 |


## Reproduce and limits

`python nesg_primary/report_split_mixed_families.py` uses existing frozen manifests and domain annotations only; it does not rerun Pfam, cluster proteins, train models or change hypotheses. **Running this script is itself an authorized-label-inspection operation**, not a normal locked-test analysis step.

Counts are descriptive. A mixed family with one example of each class is included as requested; this does not establish statistical power, assay matching, architectural identity or mechanistic evidence. The separate `overlap_sensitivity_comparison.tsv` includes exact surviving counts even when a family ceases to be mixed.
