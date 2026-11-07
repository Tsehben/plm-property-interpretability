# DeepLoc 2.0 raw data — not redistributed

These two files are **not** tracked in this repository. The DeepLoc 2.0 data page
states no license, terms of use, or redistribution grant, so we document the
download instead of redistributing it.

Download both into this directory (`DeepLoc2/raw/`):

```bash
cd DeepLoc2/raw
curl -LO https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/Swissprot_Train_Validation_dataset.csv
curl -LO https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/hpa_testset.csv
shasum -a 256 Swissprot_Train_Validation_dataset.csv hpa_testset.csv
```

Expected:

| File | Bytes | SHA-256 |
|---|---|---|
| `Swissprot_Train_Validation_dataset.csv` | 17,724,293 | `29a07b293fed2994a966b70bdcd6bacc59915b8b01fa200cb2b07d8db18384a2` |
| `hpa_testset.csv` | 570,089 | `557553bcd35e3dcb59bbfe90e1fd63341549f0a693c0e725d6c2c8387a463b37` |

Source page: <https://services.healthtech.dtu.dk/services/DeepLoc-2.0/4-Data.php>

**Reserved data.** `hpa_testset.csv` is held out and must not be read for
training or evaluation. Partition 0 of the Swiss-Prot file (2,712 records) is
likewise reserved as an internal holdout — see `../../FINDINGS.md`.
