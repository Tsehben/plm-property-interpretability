# Repo setup — read this first

MIT 6.3930 project: *Mechanistic Interpretability of Property Representations in
Protein Language Models.* Frozen ESM-2-650M, InterProt SAEs, UESolDS
(solubility) and the DeepLoc 2.0 Extracellular-vs-Cytoplasm contrast.

Start with [FINDINGS.md](FINDINGS.md) (the running lab notebook, 5 dated
entries) and [CHECKIN_1.md](CHECKIN_1.md) (the Check-in 1 brief).

---

## 1. After cloning: four things to fetch

Everything we wrote — scripts, manifests, generated cohort FASTAs, predictions,
result tables, pooled ESM-2 embeddings, reports — is **in the repo**. Four
third-party artifacts are not, for the reasons given.

| # | What | Where it goes | Why not tracked |
|---|---|---|---|
| 1 | ESM-2-650M weights (~2.6 GB) | HF cache | Model checkpoint; downloads automatically |
| 2 | InterProt SAE weights (40 MB) | `interprot_ref/raw/` | Model checkpoint; see [`interprot_ref/raw/README.md`](interprot_ref/raw/README.md) |
| 3 | DeepLoc 2.0 raw CSVs (18 MB) | `DeepLoc2/raw/` | No stated redistribution terms; see [`DeepLoc2/raw/README.md`](DeepLoc2/raw/README.md) |
| 4 | PLM_Sol trained model `model_param.t7` (28 MB) | `Violet969-PLM_Sol-3589a95/model_param/` | Upstream model weights, **unused by any of our analysis** — fetch only if you want to run PLM_Sol itself |

### 1.1 Python environment

Python 3.9 on macOS/MPS was used throughout. Exact versions recorded in the
extraction metadata:

```
transformers 4.57.6 · torch 2.8.0 · numpy 2.0.2 · scikit-learn · device=mps
```

```bash
pip install 'torch>=2.8' 'transformers>=4.57' numpy scikit-learn safetensors
```

### 1.2 ESM-2 checkpoint

Nothing to do by hand. Our scripts call
`EsmModel.from_pretrained('facebook/esm2_t33_650M_UR50D')`, which downloads to
the Hugging Face cache on first use (~2.6 GB). This is the checkpoint the
InterProt SAEs were trained on; we read `hidden_states[8]`, i.e. the output of
transformer block 8.

### 1.3 InterProt SAE

```bash
cd interprot_ref/raw
curl -L -o esm2_plm1280_l8_sae4096.safetensors \
  https://huggingface.co/liambai/InterProt-ESM2-SAEs/resolve/3ded2e51641df84d6b89f25430874e3a8eb42c24/esm2_plm1280_l8_sae4096.safetensors
shasum -a 256 esm2_plm1280_l8_sae4096.safetensors
# expect 1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd
```

Repo `liambai/InterProt-ESM2-SAEs`, revision
`3ded2e51641df84d6b89f25430874e3a8eb42c24`, Apache-2.0, 41,997,696 bytes.
The authors' reference class `sae_model.py` **is** tracked, so only the weights
are missing.

### 1.4 DeepLoc 2.0 raw data

```bash
cd DeepLoc2/raw
curl -LO https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/Swissprot_Train_Validation_dataset.csv
curl -LO https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/hpa_testset.csv
shasum -a 256 *.csv
# 29a07b293fed2994a966b70bdcd6bacc59915b8b01fa200cb2b07d8db18384a2  Swissprot_Train_Validation_dataset.csv
# 557553bcd35e3dcb59bbfe90e1fd63341549f0a693c0e725d6c2c8387a463b37  hpa_testset.csv
```

### 1.5 PLM_Sol trained model (optional)

```bash
# from https://github.com/Violet969/PLM_Sol (release V1.1)
curl -L -o /tmp/PLM_Sol-V1.1.zip https://github.com/Violet969/PLM_Sol/archive/refs/tags/V1.1.zip
unzip -j /tmp/PLM_Sol-V1.1.zip '*/model_param/model_param.t7' \
  -d Violet969-PLM_Sol-3589a95/model_param/
# expect 29,392,139 bytes
```

---

## 2. What runs without any download

The ESM-2 probes are the headline results, and both reproduce from the
tracked pooled embeddings alone — **no checkpoint, no raw CSV, no GPU**:

```bash
R="$PWD"

# UESolDS solubility probe  (~1 min)
python3 "$R/Violet969-PLM_Sol-3589a95/esm2_probe/probe_uesolds.py" \
  "$R/Violet969-PLM_Sol-3589a95" "$R/esm2_common" \
  "$R/Violet969-PLM_Sol-3589a95/esm2_probe"

# DeepLoc localization probe  (~2 min)
python3 "$R/DeepLoc2/esm2_probe/probe_deeploc.py" \
  "$R/DeepLoc2" "$R/esm2_common" "$R/DeepLoc2/esm2_probe"
```

Each rewrites its own `*_esm_results.tsv` / `.json` and OOF predictions in
place; `git diff` should come back empty.

## 3. What needs the downloads

| Step | Command | Needs |
|---|---|---|
| DeepLoc audit | `python3 DeepLoc2/audit/audit_deeploc.py DeepLoc2 DeepLoc2/audit` | 1.4 |
| DeepLoc composition baselines | `python3 DeepLoc2/baseline/baseline_localization.py DeepLoc2 DeepLoc2/baseline` | 1.4 |
| UESolDS source-vs-sequence baseline | `python3 Violet969-PLM_Sol-3589a95/baseline_source_vs_sequence/baseline.py Violet969-PLM_Sol-3589a95 Violet969-PLM_Sol-3589a95/baseline_source_vs_sequence` | tracked FASTAs only |
| UESolDS 10k/4k cohort | `python3 Violet969-PLM_Sol-3589a95/dev_cohort/build_cohort.py Violet969-PLM_Sol-3589a95 Violet969-PLM_Sol-3589a95/dev_cohort` | tracked FASTAs only |
| Rebuild `<key>\tseq` pairs | `python3 esm2_common/build_pairs.py "$PWD"` | 1.4 |
| Re-extract embeddings (~41 min + 42 min on MPS) | `python3 esm2_common/esm_extract.py --pairs <pairs.tsv> --out <emb_dir> --tag <uesolds\|deeploc> --nterm <0\|30> --dtype float16` | 1.2 |
| SAE pilot sample | `python3 esm2_common/sae_sample.py` | tracked manifests |
| SAE pilot (64+64) | `python3 esm2_common/sae_pilot.py` | 1.2 + 1.3 |

`sae_sample.py` and `sae_pilot.py` resolve the repo root as the parent of
`esm2_common/`. Override with `MIT6393_ROOT=/path/to/checkout` if needed.
Every other script takes its directories as positional arguments.

---

## 4. Reserved data — do not touch

Three holdouts are reserved and must never be read for training or evaluation.
The scripts enforce this (`sae_sample.py` asserts `partition != '0'`), and it
should stay that way.

| Reserved | Status in repo |
|---|---|
| `Violet969-PLM_Sol-3589a95/embedding_dataset/test_dataset.fasta` | tracked, unread |
| DeepLoc **partition 0** (2,712 records of the Swiss-Prot CSV) | not tracked |
| `DeepLoc2/raw/hpa_testset.csv` | not tracked |

Published FASTA/CSV originals are never modified in place.

---

## 5. Third-party licenses

| Component | License | Redistributed here? |
|---|---|---|
| PLM_Sol / UESolDS — github.com/Violet969/PLM_Sol | MIT (© 2024 Violet969), vendored at [`Violet969-PLM_Sol-3589a95/LICENSE`](Violet969-PLM_Sol-3589a95/LICENSE) | **Yes** — code + `embedding_dataset/*.fasta` |
| InterProt SAE — `liambai/InterProt-ESM2-SAEs` | Apache-2.0 | `sae_model.py` yes; 40 MB weights no (size) |
| ESM-2-650M — `facebook/esm2_t33_650M_UR50D` | MIT | No (checkpoint, auto-downloaded) |
| DeepLoc 2.0 data — DTU Health Tech | **none stated** on the data page or service page | **No** — download documented instead |

The DeepLoc call is conservative: absence of a redistribution grant is not a
grant. If someone finds an explicit license, the CSVs can be added.

---

## 6. Repo size

~206 MB tracked, dominated by the pooled ESM-2 embeddings
(`*/esm2_probe/emb/*.npz`, ~158 MB) and `Train_dataset.fasta` (22 MB). Largest
single file is 22.3 MB — comfortably under GitHub's 50 MB warning and 100 MB
hard limit, so no Git LFS is required. Clones will be slow-ish on first pull.
