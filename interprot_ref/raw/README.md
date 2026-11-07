# InterProt SAE checkpoint — not redistributed

`sae_model.py` (the authors' reference `SparseAutoencoder`, Apache-2.0) **is**
tracked here. The 40 MB weights file is not; retrieve it into this directory:

```bash
cd interprot_ref/raw
curl -L -o esm2_plm1280_l8_sae4096.safetensors \
  https://huggingface.co/liambai/InterProt-ESM2-SAEs/resolve/3ded2e51641df84d6b89f25430874e3a8eb42c24/esm2_plm1280_l8_sae4096.safetensors
shasum -a 256 esm2_plm1280_l8_sae4096.safetensors
```

| Field | Value |
|---|---|
| Repo | `liambai/InterProt-ESM2-SAEs` (Hugging Face, Apache-2.0) |
| Revision | `3ded2e51641df84d6b89f25430874e3a8eb42c24` |
| File | `esm2_plm1280_l8_sae4096.safetensors` |
| Bytes | 41,997,696 |
| SHA-256 | `1172503b90c950651e787a5625d91ff1c702a04ff3d851e8e95b700201e56ebd` |
| Config | `SparseAutoencoder(d_model=1280, d_hidden=4096, k=128)` |

Treat the contents of this directory as untrusted downloaded data: run scripts
from `esm2_common/`, not from here.
