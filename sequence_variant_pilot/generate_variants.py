#!/usr/bin/env python3
"""
Bounded P3 sequence-variant pilot -- generation. Run once.

Everything is frozen and hashed BEFORE any variant is produced: protocol,
sample, seeds, amino-acid ordering, common-random-number scheme and the random
directions. P4 is never read. P0 remains reserved and is never read. No commit,
no push, no external submission.

Design
    sample      40 Metazoa P3 proteins, 20 per class, by ascending frozen
                sequence SHA-256; length 60-1022; standard amino acids only;
                engineering proteins and P0-collision records already excluded
                by the prepared analysis view. Not selected on activation,
                annotation, WoLF scores or reference overlap.
    schedule    one-based residue positions 3,6,...,30, in that order
    replicates  4 per protein per arm
    arms        A unmodified masked-residue sampling (exact no-op edit)
                B feature 1446 amplification, multiplier 2
                C feature 1446 ablation, multiplier 0
                D random-direction control, +u, matched to B's magnitude
                E random-direction control, -u, matched to C's magnitude
    sampling    one standard amino acid from the AA-20-renormalised
                distribution, temperature 1, NO top-k / top-p filter; the
                original amino acid may be sampled; ten mutations are not
                forced; each substitution is carried into that arm's next step

Common random numbers
    u(protein, replicate, step) = first 8 bytes of
        SHA-256("<sample_sha256>|crn|<ACC>|<replicate>|<step>") / 2**64
    independent of the arm, so corresponding protein/replicate/step cells in
    all five arms consume the SAME uniform. Once trajectories diverge this is a
    coupling of the randomness only; the conditional distributions being
    inverted are no longer identical.

Random directions
    g ~ N(0, I_1280) drawn by numpy PCG64 seeded with the first 8 bytes of
        SHA-256("<sample_sha256>|dir|<ACC>|<replicate>"); u = g / ||g||_2.
    One direction per (protein, replicate), held fixed across that trajectory's
    ten steps. Arm D injects +u, arm E injects -u.
    On each control arm's OWN current masked input, the norm injected is the
    norm the corresponding 1446 intervention WOULD produce on THAT input. If
    1446 is inactive there, the perturbation is exactly zero. The injection is a
    raw residual-space vector: SAE top-K support is never modified to implement
    it.

Usage:  python3 generate_variants.py
"""

import hashlib
import json
import os
import platform
import sys
import time
from collections import Counter

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "v1")
PER = os.path.join(OUT, "per_protein")
sys.path.insert(0, os.path.join(ROOT, "interprot_ref", "raw"))
sys.path.insert(0, os.path.join(ROOT, "candidate_discovery"))
sys.path.insert(0, os.path.join(ROOT, "intervention_pilot"))
from sae_model import SparseAutoencoder  # noqa: E402
import run_discovery as RD  # noqa: E402
import run_pilot as RP  # noqa: E402

REF_SHA_INDEX = os.path.join(ROOT, "wolfpsort_setup", "reference_audit",
                             "reference_sequence_sha256.tsv")

N_PER_CLASS = 20
KINGDOM = "Metazoa"
POSITIONS = [3, 6, 9, 12, 15, 18, 21, 24, 27, 30]
N_REPLICATES = 4
CANDIDATE = 1446
AA_ORDER = sorted(RP.STANDARD_AA)          # canonical ordering: A C D E F G H I K L M N P Q R S T V W Y
TEMPERATURE = 1.0

# arm -> (kind, multiplier, direction sign)
ARMS = {
    "A": ("unmodified", 1.0, 0),
    "B": ("candidate", 2.0, 0),
    "C": ("candidate", 0.0, 0),
    "D": ("random", 2.0, +1),
    "E": ("random", 0.0, -1),
}

CONFIG = dict(
    purpose="bounded exploratory P3 sequence-variant pilot, run once at the frozen "
            "sample and replicate counts",
    partition="P3 (tuning) only. P4 is not read. P0 remains reserved and is never read.",
    candidate=CANDIDATE,
    candidate_provenance=(
        "feature 1446 was selected during discovery on P1+P2 as a representative "
        "signal-associated feature (rank 3, AUROC 0.974, higher in Extracellular), before "
        "any P4 output existed. No new candidate selection is performed here, and the P4 "
        "result for 2564 is not used to choose a candidate."
    ),
    n_per_class=N_PER_CLASS, kingdom=KINGDOM,
    selection_rule=(
        "within each localisation class, sort eligible P3 Metazoa proteins by frozen "
        "sequence SHA-256 hex ascending and take the first 20. No random seed. No "
        "filtering on feature activation, annotation, WoLF PSORT scores or reference "
        "overlap."
    ),
    eligibility="role == tuning and in_core == 1 in the prepared inclusion manifest "
                "(engineering proteins and p0_base_collision records already excluded); "
                "kingdom_deeploc == Metazoa; original length 60-1022; standard 20 amino "
                "acids only",
    kingdom_rationale=(
        "one WoLF PSORT organism model (animal) must apply to every protein: k, the "
        "reference set and the class vocabulary all differ between models, so scores are "
        "not comparable across them and per-protein model switching is not an option"
    ),
    positions=POSITIONS, positions_note="one-based residue coordinates, applied in this "
                                        "order, one masked position per step; position 1 "
                                        "is never resampled",
    n_replicates=N_REPLICATES,
    arms={k: dict(kind=v[0], multiplier=v[1], direction_sign=v[2]) for k, v in ARMS.items()},
    edit="validated direct decoder-delta at the currently masked token only, hook on "
         "esm.encoder.layer[7]",
    sampling=dict(
        distribution="AA-20 renormalised softmax at the masked token of the arm's own "
                     "current forward pass",
        temperature=TEMPERATURE,
        filters="none: no top-k, no top-p, no repetition penalty, no banned tokens",
        amino_acid_ordering=AA_ORDER,
        method="inverse CDF over AA_ORDER using the shared uniform",
        original_allowed="yes; sampling the original amino acid is permitted and is not "
                         "resampled. Ten mutations are NOT forced.",
        carry_forward="each sampled substitution is written into that arm's sequence and "
                      "carried into its next step",
    ),
    common_random_numbers=dict(
        scheme="u = int(SHA-256('<sample_sha256>|crn|<ACC>|<replicate>|<step>')[:8]) / 2**64",
        arm_independent=True,
        caveat="identical uniforms across arms. After the trajectories diverge this is a "
               "coupling of the randomness only -- the conditional distributions being "
               "inverted are no longer identical, so this is variance reduction, not exact "
               "pairing.",
    ),
    random_directions=dict(
        distribution="isotropic on the unit sphere S^1279: g ~ N(0, I_1280), u = g/||g||_2",
        generator="numpy PCG64 (np.random.default_rng) seeded per (protein, replicate)",
        seed="int(SHA-256('<sample_sha256>|dir|<ACC>|<replicate>')[:8])",
        n_directions="one per (protein, replicate) = 4 per protein, held fixed across that "
                     "trajectory's ten steps",
        signs="arm D injects +u, arm E injects -u",
        norm_rule="on each control arm's OWN current masked input, inject the norm that the "
                  "corresponding 1446 intervention would produce on THAT input; if 1446 is "
                  "inactive there, inject exactly zero",
        matching_description="a matched intervention rule on the current input. It is NOT "
                             "perfect matching of entire sequence-generation trajectories: "
                             "after arms diverge, paired arms need not have identical norms, "
                             "and both the hypothetical control norms and the realised "
                             "candidate-arm norms are recorded separately.",
        top_k_policy="the random direction is injected as a raw residual-space vector. SAE "
                     "top-K support is never modified to implement it.",
        control_caveat=(
            "these are NOT the activation-matched SAE controls. Under the prespecified "
            "+/-20% rule feature 1446 had 0 of 5 eligible matched SAE controls. A random "
            "isotropic direction matches no activation statistic and lies outside the "
            "sparse code; it controls for perturbation magnitude at the edited token only."
        ),
    ),
    readout=dict(
        name="margin difference D",
        margin="margin(s) = single-site score(extr) - single-site score(cyto)",
        D="D = margin(variant) - margin(original)",
        organism_model="animal, fixed for every protein",
        extraction="PENDING VERIFICATION in score_variants.py: the summary line omits "
                   "nonzero classes, so scores are recomputed from the --print-neighbors "
                   "list through data/animal.wolfu and asserted to reproduce every printed "
                   "score exactly",
        naming="predictor scores. Not probabilities and not biological measurements.",
    ),
    compound_class_rule=(
        "OPERATIONALLY DEFINED HERE, BEFORE GENERATION: compound (dual-site) classes are "
        "never added into the extr or cyto terms, because score('extr') already carries "
        "half-credit for every extr_* neighbour and adding a compound score would double "
        "count the same neighbours. Compound outputs are reported descriptively via the "
        "flags top_is_compound / top_mentions_extr / top_mentions_cyto. No compound-class "
        "sensitivity rule beyond this is defined, so none will be introduced after seeing "
        "results."
    ),
    analysis=dict(
        unit="protein. Replicates are averaged within protein before any across-protein "
             "summary. Steps, arms and replicates are not independent sample units.",
        contrasts=["each arm's paired original-to-variant D", "B-A", "C-A", "B-D", "C-E"],
        status="exploratory P3 pilot evidence. No efficacy pass threshold, no confirmatory "
               "claim, no decision taken on whether any interval contains zero.",
    ),
    prohibited="no regeneration of unfavourable variants; no P4 or P0 access; no new "
               "candidate selection; no external sequence submission; no commit or push",
)


def config_hash():
    return hashlib.sha256(json.dumps(CONFIG, sort_keys=True, default=str).encode()).hexdigest()


def code_version():
    h = {}
    for p in [os.path.abspath(__file__),
              os.path.join(ROOT, "intervention_pilot", "run_pilot.py"),
              os.path.join(ROOT, "candidate_discovery", "run_discovery.py"),
              os.path.join(ROOT, "interprot_ref", "raw", "sae_model.py")]:
        h[os.path.relpath(p, ROOT)] = RD.sha256_file(p)
    return h


def crn_uniform(sample_hash, acc, rep, step):
    d = hashlib.sha256(f"{sample_hash}|crn|{acc}|{rep}|{step}".encode()).digest()
    return int.from_bytes(d[:8], "big") / float(1 << 64)


def direction_seed(sample_hash, acc, rep):
    d = hashlib.sha256(f"{sample_hash}|dir|{acc}|{rep}".encode()).digest()
    return int.from_bytes(d[:8], "big")


def unit_direction(seed, dim=RD.HIDDEN):
    g = np.random.default_rng(seed).standard_normal(dim)
    return g / np.linalg.norm(g)


def kd_mean(seq):
    return float(np.mean([RP.KD[a] for a in seq]))


def build_sample():
    import csv
    incl = {r["ACC"]: r for r in RD.read_tsv(RP.INCLUSION)}
    manifest = {r["ACC"]: r for r in RD.read_tsv(RP.MANIFEST)}
    seqs = {}
    with open(RP.RAW_CSV, newline="") as fh:
        for r in csv.DictReader(fh):
            seqs[r["ACC"]] = r["Sequence"]
    pool, excl = [], Counter()
    for acc, r in incl.items():
        if r["role"] != "tuning":
            excl["not_P3"] += 1
            continue
        if r["in_core"] != "1":
            excl[f"excluded_by_view:{r['exclusion_reason']}"] += 1
            continue
        m = manifest[acc]
        assert m["partition"] == "3", acc
        if m["kingdom"] != KINGDOM:
            excl[f"kingdom_{m['kingdom']}"] += 1
            continue
        s = seqs[acc]
        assert hashlib.sha256(s.encode()).hexdigest() == m["sha256"], acc
        if not (RP.LEN_MIN <= len(s) <= RP.LEN_MAX):
            excl["length_outside_60_1022"] += 1
            continue
        if set(s) - RP.STANDARD_AA:
            excl["nonstandard_residues"] += 1
            continue
        pool.append(dict(ACC=acc, seq=s, target=int(m["target"]), original_len=len(s),
                         kingdom=m["kingdom"], sha256=m["sha256"]))
    chosen, quota = [], {}
    for t in (1, 0):
        cls = sorted([p for p in pool if p["target"] == t], key=lambda p: p["sha256"])
        quota[str(t)] = dict(available=len(cls), requested=N_PER_CLASS)
        chosen += cls[:N_PER_CLASS]
    for t, q in quota.items():
        assert q["available"] >= N_PER_CLASS, (
            f"QUOTA CANNOT BE MET for class {t}: {q}. Stopping rather than changing the rule.")
    assert len(chosen) == 2 * N_PER_CLASS
    assert len({p["ACC"] for p in chosen}) == len(chosen)
    return chosen, dict(pool_size=len(pool), quota=quota, exclusions=dict(excl))


def reference_overlap(sample):
    """Exact-sequence overlap against the shipped WoLF PSORT reference sequences.
    Descriptive only; the frozen sample is retained whatever this returns."""
    idx = {}
    with open(REF_SHA_INDEX) as fh:
        head = fh.readline()
        assert head.split("\t")[0] == "sha256"
        for ln in fh:
            h, org, ident, cls, length = ln.rstrip("\n").split("\t")
            idx.setdefault(h, []).append(dict(organism=org, reference_id=ident,
                                              reference_class=cls, length=int(length)))
    rows = []
    for p in sample:
        h = hashlib.sha256(p["seq"].encode()).hexdigest()
        hits = idx.get(h, [])
        rows.append(dict(
            ACC=p["ACC"], localisation_class="Extracellular" if p["target"] else "Cytoplasm",
            original_len=p["original_len"], sequence_sha256=p["sha256"],
            exact_reference_match=int(bool(hits)), n_matches=len(hits),
            reference_ids=";".join(x["reference_id"] for x in hits),
            reference_organism_models=";".join(sorted({x["organism"] for x in hits})),
            reference_classes=";".join(x["reference_class"] for x in hits)))
    return rows


class MulEdit:
    """Direct decoder-delta edit of one latent coordinate at one token, specified
    by a multiplier on the feature's own clean activation. multiplier == 1.0 is an
    exact no-op (a_ed is bit-identical to a, so the delta is exactly zero)."""

    def __init__(self, sae, feature, pos, multiplier):
        self.sae, self.feature, self.pos, self.mult = sae, feature, pos, multiplier
        self.info = {}

    def __call__(self, module, inp, out):
        is_tuple = isinstance(out, tuple)
        h = out[0] if is_tuple else out
        x = h[0].float().cpu()
        sae = self.sae
        a = sae.get_acts(x)
        before = float(a[self.pos, self.feature])
        new = before * self.mult
        a_ed = a.clone()
        a_ed[self.pos, self.feature] = new
        std = sae.encode(x)[2]
        ta, te = sae.topK_activation(a, sae.k), sae.topK_activation(a_ed, sae.k)
        delta = ((te - ta) @ sae.w_dec) * std
        other = torch.ones(sae.d_hidden, dtype=torch.bool)
        other[self.feature] = False
        sup_a, sup_e = ta > 0, te > 0
        self.info = dict(
            arm_kind="unmodified" if self.mult == 1.0 else "candidate",
            candidate_latent_clean=before, latent_after=new, latent_delta=new - before,
            n_latent_coords_changed=int((a_ed != a).sum()),
            support_evicted_unintended=int((sup_a & ~sup_e & other).sum()),
            support_inserted_unintended=int((~sup_a & sup_e & other).sum()),
            delta_h_l2=float(delta.norm()),
            n_nonzero_delta_positions=int((delta.abs().sum(-1) > 0).sum()),
            delta_confined_to_edited_position=int(
                torch.unique((delta.abs().sum(-1) > 0).nonzero().flatten()).tolist()
                in ([], [self.pos])),
            std_at_pos=float(std[self.pos]),
            h_l2_norm=float(x.norm()),
        )
        new_h = h.clone()
        new_h[0] = (x + delta).to(h.dtype).to(h.device)
        return (new_h,) + out[1:] if is_tuple else new_h


class RandEdit:
    """Inject a fixed unit residual-space direction at one token, scaled to the
    norm the 1446 intervention WOULD produce on THIS input. The SAE latents and
    the top-K support are read but never modified."""

    def __init__(self, sae, feature, pos, multiplier, direction):
        self.sae, self.feature, self.pos, self.mult = sae, feature, pos, multiplier
        self.direction = torch.as_tensor(direction, dtype=torch.float32)
        self.info = {}

    def __call__(self, module, inp, out):
        is_tuple = isinstance(out, tuple)
        h = out[0] if is_tuple else out
        x = h[0].float().cpu()
        sae = self.sae
        a = sae.get_acts(x)
        before = float(a[self.pos, self.feature])
        a_ed = a.clone()
        a_ed[self.pos, self.feature] = before * self.mult
        std = sae.encode(x)[2]
        cand_delta = ((sae.topK_activation(a_ed, sae.k)
                       - sae.topK_activation(a, sae.k)) @ sae.w_dec) * std
        T = float(cand_delta.norm())
        delta = torch.zeros_like(x)
        delta[self.pos] = self.direction * T
        realised = float(delta.norm())
        self.info = dict(
            arm_kind="random",
            candidate_latent_clean=before,
            hypothetical_candidate_norm_on_this_input=T,
            injected_norm=realised,
            rel_norm_error=(abs(realised - T) / T) if T > 0 else 0.0,
            cosine_with_candidate_delta=(
                float(torch.dot(delta[self.pos], cand_delta[self.pos])
                      / (realised * T)) if T > 0 and realised > 0 else 0.0),
            n_latent_coords_changed=0,
            support_evicted_unintended=0, support_inserted_unintended=0,
            delta_h_l2=realised,
            n_nonzero_delta_positions=int((delta.abs().sum(-1) > 0).sum()),
            delta_confined_to_edited_position=int(
                torch.unique((delta.abs().sum(-1) > 0).nonzero().flatten()).tolist()
                in ([], [self.pos])),
            std_at_pos=float(std[self.pos]),
            h_l2_norm=float(x.norm()),
        )
        new_h = h.clone()
        new_h[0] = (x + delta).to(h.dtype).to(h.device)
        return (new_h,) + out[1:] if is_tuple else new_h


def main():
    t0 = time.time()
    os.makedirs(PER, exist_ok=True)

    sample, pool_info = build_sample()
    sample_hash = hashlib.sha256(
        "\n".join(f"{p['ACC']}\t{p['sha256']}\t{p['target']}" for p in sample).encode()
    ).hexdigest()

    dirs = {(p["ACC"], r): direction_seed(sample_hash, p["ACC"], r)
            for p in sample for r in range(N_REPLICATES)}
    overlap = reference_overlap(sample)

    frozen = dict(
        config=CONFIG, config_hash=config_hash(), code_version=code_version(),
        sample=dict(n=len(sample), sample_sha256=sample_hash,
                    n_Extracellular=sum(1 for p in sample if p["target"] == 1),
                    n_Cytoplasm=sum(1 for p in sample if p["target"] == 0),
                    accessions=[p["ACC"] for p in sample],
                    sequence_sha256={p["ACC"]: p["sha256"] for p in sample},
                    **pool_info),
        direction_seeds={f"{a}|{r}": s for (a, r), s in sorted(dirs.items())},
        reference_overlap_summary=dict(
            n_exact_matches=sum(r["exact_reference_match"] for r in overlap),
            note="descriptive; the frozen sample is retained whatever this shows. Exact "
                 "nonmatches are NOT evidence of homology independence."),
        inputs=dict(inclusion_manifest=RD.sha256_file(RP.INCLUSION),
                    deeploc_manifest=RD.sha256_file(RP.MANIFEST),
                    selected_feature_ids=RD.sha256_file(RP.SELECTED),
                    reference_sequence_index=RD.sha256_file(REF_SHA_INDEX),
                    sae_file=RD.SAE_SHA),
        environment=dict(python=platform.python_version(), torch=torch.__version__,
                         transformers=__import__("transformers").__version__,
                         numpy=np.__version__, device=str(RP.device_of()),
                         platform=platform.platform()),
        frozen_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    )
    fc = os.path.join(OUT, "frozen_protocol.json")
    if os.path.exists(fc):
        old = json.load(open(fc))
        assert old["sample"]["sample_sha256"] == sample_hash, "sample drifted after freezing"
        assert old["config_hash"] == config_hash(), "protocol drifted after freezing"
        print("frozen protocol already present and matches", flush=True)
    else:
        with open(fc, "w") as fh:
            json.dump(frozen, fh, indent=2)
    RD.write_tsv(os.path.join(OUT, "pilot_sample.tsv"),
                 ["ACC", "target", "localisation_class", "original_len", "kingdom", "sha256"],
                 [dict(ACC=p["ACC"], target=p["target"],
                       localisation_class="Extracellular" if p["target"] else "Cytoplasm",
                       original_len=p["original_len"], kingdom=p["kingdom"], sha256=p["sha256"])
                  for p in sample])
    RD.write_tsv(os.path.join(OUT, "reference_overlap.tsv"), list(overlap[0].keys()), overlap)
    print(f"FROZEN  protocol={config_hash()[:16]}  sample={sample_hash[:16]}  "
          f"n={len(sample)} ({frozen['sample']['n_Extracellular']}/"
          f"{frozen['sample']['n_Cytoplasm']})  "
          f"exact_reference_matches={frozen['reference_overlap_summary']['n_exact_matches']}",
          flush=True)

    # ------------------------------------------------------------------ models
    dev = RP.device_of()
    tok = AutoTokenizer.from_pretrained(RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION)
    mlm = EsmForMaskedLM.from_pretrained(RD.CHECKPOINT, revision=RD.CHECKPOINT_REVISION,
                                         torch_dtype=torch.float32)
    mlm.eval().to(dev)
    for prm in mlm.parameters():
        prm.requires_grad_(False)
    from safetensors.torch import load_file
    assert RD.sha256_file(RD.SAE_FILE) == RD.SAE_SHA
    sae = SparseAutoencoder(RD.HIDDEN, RD.SAE_DIM)
    sae.load_state_dict(load_file(RD.SAE_FILE))
    sae.eval()
    aa_ids = [tok.convert_tokens_to_ids(a) for a in AA_ORDER]
    assert len(set(aa_ids)) == 20
    mask_id = tok.mask_token_id

    n_pass = 0
    for pi, prot in enumerate(sample):
        fp = os.path.join(PER, f"{prot['ACC']}.json")
        if os.path.exists(fp):
            continue
        acc, orig, L = prot["ACC"], prot["seq"], prot["original_len"]
        steps, variants = [], []
        for rep in range(N_REPLICATES):
            u_vec = unit_direction(dirs[(acc, rep)])
            for arm, (kind, mult, sign) in ARMS.items():
                seq = orig
                for si, p1 in enumerate(POSITIONS):
                    enc = tok([seq], return_tensors="pt", add_special_tokens=True)
                    ids0, am = enc["input_ids"], enc["attention_mask"]
                    assert int(am.sum()) == L + 2
                    tpos = p1
                    tokname = tok.convert_ids_to_tokens(int(ids0[0, tpos]))
                    assert tokname == seq[p1 - 1], (acc, arm, rep, si, tokname, seq[p1 - 1])
                    ids = ids0.clone()
                    ids[0, tpos] = mask_id
                    idsd, amd = ids.to(dev), am.to(dev)
                    if kind == "random":
                        hk = RandEdit(sae, CANDIDATE, tpos, mult, u_vec * sign)
                    else:
                        hk = MulEdit(sae, CANDIDATE, tpos, mult)
                    logits = RP.forward_with(mlm, idsd, amd, hk)
                    n_pass += 1
                    q = torch.softmax(
                        logits[0, tpos, aa_ids].cpu().double() / TEMPERATURE, -1).numpy()
                    u = crn_uniform(sample_hash, acc, rep, si)
                    cdf = np.cumsum(q)
                    j = int(np.searchsorted(cdf, u, side="right"))
                    j = min(j, 19)
                    aa_new = AA_ORDER[j]
                    aa_old = seq[p1 - 1]
                    steps.append(dict(
                        ACC=acc, target=prot["target"], arm=arm,
                        replicate=rep, step=si, residue_position=p1, token_index=tpos,
                        multiplier=mult, direction_sign=sign,
                        residue_before=aa_old, residue_sampled=aa_new,
                        changed=int(aa_new != aa_old),
                        uniform=f"{u:.17g}",
                        prob_sampled=f"{float(q[j]):.17g}",
                        prob_original=f"{float(q[AA_ORDER.index(aa_old)]):.17g}",
                        top1_aa=AA_ORDER[int(q.argmax())],
                        prob_top1=f"{float(q.max()):.17g}",
                        entropy_nats=f"{float(-(q * np.log(q + 1e-300)).sum()):.17g}",
                        **{k: (f"{v:.17g}" if isinstance(v, float) else v)
                           for k, v in hk.info.items()}))
                    seq = seq[:p1 - 1] + aa_new + seq[p1:]
                assert len(seq) == L
                assert not (set(seq) - RP.STANDARD_AA)
                diffs = [i for i in range(L) if seq[i] != orig[i]]
                assert all((i + 1) in POSITIONS for i in diffs), (acc, arm, rep, diffs)
                variants.append(dict(
                    ACC=acc, target=prot["target"], arm=arm, replicate=rep,
                    variant_id=f"{acc}__{arm}__r{rep}", sequence=seq,
                    sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                    length=len(seq), n_mutations=len(diffs),
                    mutated_positions=";".join(str(i + 1) for i in diffs),
                    mutations=";".join(f"{orig[i]}{i+1}{seq[i]}" for i in diffs),
                    identical_to_original=int(seq == orig),
                    sequence_identity=f"{(L - len(diffs)) / L:.17g}",
                    kd_whole=f"{kd_mean(seq):.17g}", kd_nterm30=f"{kd_mean(seq[:30]):.17g}",
                    delta_kd_whole=f"{kd_mean(seq) - kd_mean(orig):.17g}",
                    delta_kd_nterm30=f"{kd_mean(seq[:30]) - kd_mean(orig[:30]):.17g}",
                    direction_seed=dirs[(acc, rep)]))
        with open(fp, "w") as fh:
            json.dump(dict(protein=acc, target=prot["target"], original_len=L,
                           original_sequence=orig,
                           original_sha256=prot["sha256"],
                           kd_whole_original=kd_mean(orig),
                           kd_nterm30_original=kd_mean(orig[:30]),
                           positions=POSITIONS, variants=variants, steps=steps), fh)
        el = time.time() - t0
        print(f"  [{pi+1}/{len(sample)}] {acc} passes={n_pass} {el:.0f}s "
              f"({n_pass/max(el,1e-9):.2f} pass/s)", flush=True)

    with open(os.path.join(OUT, "generation_report.json"), "w") as fh:
        json.dump(dict(forward_passes=n_pass, seconds=time.time() - t0,
                       n_proteins=len(sample), n_replicates=N_REPLICATES,
                       n_arms=len(ARMS), n_steps=len(POSITIONS),
                       expected_passes=len(sample) * N_REPLICATES * len(ARMS) * len(POSITIONS)),
                  fh, indent=2)
    print(f"DONE  {n_pass} passes  {time.time()-t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
