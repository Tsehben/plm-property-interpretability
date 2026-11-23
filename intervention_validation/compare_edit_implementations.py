#!/usr/bin/env python3
"""
Compare two implementations of the error-preserving SAE feature edit, against a
float64 reference, on the exact inputs/features/strengths of the completed
validation.

  DoR  difference-of-reconstructions   delta = decode(a_edited, mu, std) - decode(a, mu, std)
  DD   direct decoder-delta            delta = ((topK(a_edited) - topK(a)) @ w_dec) * std

Both are algebraically identical: b_pre and mu cancel in the DoR difference.
They differ only in floating-point behaviour. DoR forms two large tensors and
subtracts them; DD forms the (sparse) latent difference first.

Writes out/edit_implementation_comparison.json. Does not modify any existing
output. Does not change edit strengths.

Run:  python3 intervention_validation/compare_edit_implementations.py
"""
import csv
import json
import os
import sys
import time

import numpy as np
import torch
from transformers import AutoTokenizer, EsmForMaskedLM

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF = os.path.join(ROOT, 'interprot_ref', 'raw')
OUT = os.path.join(ROOT, 'intervention_validation', 'out')
sys.path.insert(0, REF)
from sae_model import SparseAutoencoder  # noqa: E402

CHECKPOINT = 'facebook/esm2_t33_650M_UR50D'
SAE_FILE = os.path.join(REF, 'esm2_plm1280_l8_sae4096.safetensors')
HIDDEN, SAE_DIM, TOPK = 1280, 4096, 128
MAX_RESIDUES, N_SAMPLE = 1022, 8
HOOK_LAYER = 7                 # verified in validate_intervention.py
AMPLIFY_FACTOR = 2.0           # unchanged from the validation

# Elementwise relative error is only meaningful where the reference is not
# itself ~0. TAU_REL defines "not ~0" as a fraction of the reference peak.
TAU_REL = 1e-3


def delta_DoR(sae, x, a, a_ed):
    """Difference of two full reconstructions."""
    _, mu, std = sae.encode(x)
    return sae.decode(a_ed, mu, std) - sae.decode(a, mu, std)


def delta_DD(sae, x, a, a_ed):
    """Direct decoder-delta: ((topK(a_ed) - topK(a)) @ w_dec) * std."""
    _, _, std = sae.encode(x)
    return ((sae.topK_activation(a_ed, sae.k)
             - sae.topK_activation(a, sae.k)) @ sae.w_dec) * std


def err_metrics(approx, ref):
    """Absolute and explicitly-defined relative errors. ref is float64."""
    e = (approx.double() - ref)
    abs_max = float(e.abs().max())
    abs_l2 = float(e.norm())
    ref_l2 = float(ref.norm())
    ref_peak = float(ref.abs().max())

    out = {
        'abs_err_max': abs_max,
        'abs_err_l2': abs_l2,
        'ref_l2_norm': ref_l2,
        'ref_peak_abs': ref_peak,
        # Primary relative measure: no denominator pathology unless the whole
        # reference delta is identically zero.
        'rel_err_l2': (abs_l2 / ref_l2) if ref_l2 > 0 else None,
        'rel_err_max_vs_peak': (abs_max / ref_peak) if ref_peak > 0 else None,
    }

    # Elementwise relative error, restricted to coordinates where the reference
    # magnitude is at least TAU_REL * peak. Coordinates below that threshold are
    # excluded and counted, rather than producing an inflated ratio.
    if ref_peak > 0:
        tau = TAU_REL * ref_peak
        keep = ref.abs() >= tau
        n_keep = int(keep.sum())
        if n_keep:
            rel = (e[keep].abs() / ref[keep].abs())
            out.update({
                'elementwise_tau_abs': tau,
                'elementwise_n_included': n_keep,
                'elementwise_n_excluded': int(ref.numel() - n_keep),
                'elementwise_rel_err_max': float(rel.max()),
                'elementwise_rel_err_median': float(rel.median()),
            })
    # Digits of agreement implied by the RELATIVE error (never from abs alone).
    r = out['rel_err_l2']
    out['significant_digits_agreeing_l2'] = (float(-np.log10(r)) if r and r > 0
                                             else ('exact' if r == 0 else None))
    return out


def main():
    t0 = time.time()
    dev = torch.device('mps') if torch.backends.mps.is_available() else torch.device('cpu')
    from safetensors.torch import load_file
    sd = load_file(SAE_FILE)

    sae32 = SparseAutoencoder(HIDDEN, SAE_DIM)
    sae32.load_state_dict(sd)
    sae32.eval()
    sae64 = SparseAutoencoder(HIDDEN, SAE_DIM)
    sae64.load_state_dict(sd)
    sae64 = sae64.double()
    sae64.eval()
    for m in (sae32, sae64):
        for p in m.parameters():
            p.requires_grad_(False)

    tok = AutoTokenizer.from_pretrained(CHECKPOINT)
    mlm = EsmForMaskedLM.from_pretrained(CHECKPOINT, dtype=torch.float32)
    mlm.eval().to(dev)
    for p in mlm.parameters():
        p.requires_grad_(False)

    # --- same engineering sample, same rule as the validation ---------------
    with open(os.path.join(ROOT, 'DeepLoc2', 'sae_pilot', 'sae_pilot_sample.tsv'),
              newline='') as f:
        rows = [r for r in csv.DictReader(f, delimiter='\t')
                if int(r['seq_len']) <= MAX_RESIDUES]
    rows.sort(key=lambda r: r['key'])
    sample = [(r['key'], r['sequence']) for r in rows[:N_SAMPLE]]

    results = []
    support_mismatch = 0

    for key, seq in sample:
        L = len(seq)
        mask_tok, editB_tok = L // 2 + 1, L // 4 + 1
        enc = tok([seq], return_tensors='pt')
        ids = enc['input_ids'].to(dev).clone()
        am = enc['attention_mask'].to(dev)
        ids[0, mask_tok] = tok.mask_token_id

        cap = {}

        def grab(mod, inp, out):
            t = out[0] if isinstance(out, tuple) else out
            cap['h'] = t.detach().clone()

        hd = mlm.esm.encoder.layer[HOOK_LAYER].register_forward_hook(grab)
        with torch.inference_mode():
            clean_logits = mlm(input_ids=ids, attention_mask=am).logits.detach().clone()
        hd.remove()

        h = cap['h']
        x32 = h[0].float().cpu()
        x64 = x32.double()
        a32 = sae32.get_acts(x32)
        a64 = sae64.get_acts(x64)

        # Does top-k select the same support in fp32 and fp64?
        if not bool(((a32 > 0) == (a64 > 0)).all()):
            support_mismatch += 1

        for case, edit_tok in [('A_edit_at_masked_position', mask_tok),
                               ('B_edit_at_distinct_unmasked_position', editB_tok)]:
            act = a32[edit_tok]
            live = (act > 0).nonzero().flatten()
            order = sorted(live.tolist(), key=lambda j: (-float(act[j]), j))
            rules = {'strongest_active': order[0],
                     'median_rank_active': order[len(order) // 2]}

            for rule, feat in rules.items():
                for arm, scale in [('zero_strength_edit', 1.0),
                                   ('amplify_x2', AMPLIFY_FACTOR),
                                   ('ablate_to_zero', 0.0)]:
                    a32e = a32.clone()
                    a32e[edit_tok, feat] = a32[edit_tok, feat] * scale
                    a64e = a64.clone()
                    a64e[edit_tok, feat] = a64[edit_tok, feat] * scale

                    d_dor32 = delta_DoR(sae32, x32, a32, a32e)
                    d_dd32 = delta_DD(sae32, x32, a32, a32e)
                    d_dor64 = delta_DoR(sae64, x64, a64, a64e)
                    d_dd64 = delta_DD(sae64, x64, a64, a64e)

                    # Reference: float64 direct form. Justified by showing the
                    # two float64 forms agree far below fp32 resolution.
                    ref = d_dd64
                    ref_l2 = float(ref.norm())
                    f64_cross = {
                        'abs_max': float((d_dor64 - d_dd64).abs().max()),
                        'rel_l2': (float((d_dor64 - d_dd64).norm() / ref_l2)
                                   if ref_l2 > 0 else None),
                    }

                    rec = {
                        'ACC': key, 'case': case, 'rule': rule, 'arm': arm,
                        'feature': feat, 'scale': scale,
                        'edit_token_index': edit_tok,
                        'ref_delta_l2': ref_l2,
                        'float64_DoR_vs_DD': f64_cross,
                        'DoR_fp32_vs_ref64': err_metrics(d_dor32, ref),
                        'DD_fp32_vs_ref64': err_metrics(d_dd32, ref),
                        'DoR_vs_DD_fp32': {
                            'abs_max': float((d_dor32 - d_dd32).abs().max()),
                            'bitwise_equal': bool(torch.equal(d_dor32, d_dd32)),
                        },
                    }

                    # No-op arms: exactness is the whole point.
                    if scale == 1.0:
                        rec['noop'] = {
                            'DoR_delta_exactly_zero': bool((d_dor32 == 0).all()),
                            'DD_delta_exactly_zero': bool((d_dd32 == 0).all()),
                            'DoR_delta_max_abs': float(d_dor32.abs().max()),
                            'DD_delta_max_abs': float(d_dd32.abs().max()),
                        }

                    # Delta localisation: nonzero rows must be the edit position.
                    for nm, d in (('DoR', d_dor32), ('DD', d_dd32)):
                        nz = torch.unique(
                            (d.abs().sum(-1) > 0).nonzero().flatten()).tolist()
                        rec[f'{nm}_nonzero_rows'] = nz
                        rec[f'{nm}_delta_confined'] = nz in ([], [edit_tok])

                    # Downstream: inject each fp32 delta and compare logits.
                    logits = {}
                    for nm, d in (('DoR', d_dor32), ('DD', d_dd32),
                                  ('ref64_rounded', (x64 + ref).float() - x32)):
                        def mk(dd):
                            def hook(mod, inp, out):
                                is_t = isinstance(out, tuple)
                                hh = out[0] if is_t else out
                                nh = hh.clone()
                                nh[0] = (hh[0].float().cpu() + dd).to(hh.dtype).to(hh.device)
                                return (nh,) + out[1:] if is_t else nh
                            return hook
                        hdl = mlm.esm.encoder.layer[HOOK_LAYER].register_forward_hook(mk(d))
                        with torch.inference_mode():
                            logits[nm] = mlm(input_ids=ids,
                                             attention_mask=am).logits.detach().clone()
                        hdl.remove()

                    def ldiff(p, q):
                        return {'max_abs_logit_diff': float((p - q).abs().max()),
                                'bitwise_equal': bool(torch.equal(p, q)),
                                'max_abs_at_mask': float(
                                    (p[0, mask_tok] - q[0, mask_tok]).abs().max())}

                    rec['downstream'] = {
                        'DoR_vs_DD': ldiff(logits['DoR'], logits['DD']),
                        'DoR_vs_ref64rounded': ldiff(logits['DoR'], logits['ref64_rounded']),
                        'DD_vs_ref64rounded': ldiff(logits['DD'], logits['ref64_rounded']),
                        'DoR_vs_clean': ldiff(logits['DoR'], clean_logits),
                        'DD_vs_clean': ldiff(logits['DD'], clean_logits),
                    }
                    results.append(rec)
        print(f'  {key}: done ({len(results)} conditions so far)', flush=True)

    # ---------------- aggregate ----------------
    nz = [r for r in results if r['scale'] != 1.0]
    noop = [r for r in results if r['scale'] == 1.0]

    def col(rs, path):
        out = []
        for r in rs:
            v = r
            for k in path:
                v = v[k]
            if v is not None:
                out.append(v)
        return out

    agg = {
        'n_conditions_total': len(results),
        'n_nonzero_edit_conditions': len(nz),
        'n_noop_conditions': len(noop),
        'fp32_fp64_topk_support_mismatches': support_mismatch,
        'float64_DoR_vs_DD_rel_l2_max': max(col(nz, ['float64_DoR_vs_DD', 'rel_l2'])),
        'DoR_fp32_rel_err_l2': {
            'max': max(col(nz, ['DoR_fp32_vs_ref64', 'rel_err_l2'])),
            'median': float(np.median(col(nz, ['DoR_fp32_vs_ref64', 'rel_err_l2']))),
        },
        'DD_fp32_rel_err_l2': {
            'max': max(col(nz, ['DD_fp32_vs_ref64', 'rel_err_l2'])),
            'median': float(np.median(col(nz, ['DD_fp32_vs_ref64', 'rel_err_l2']))),
        },
        'DoR_fp32_abs_err_max': max(col(nz, ['DoR_fp32_vs_ref64', 'abs_err_max'])),
        'DD_fp32_abs_err_max': max(col(nz, ['DD_fp32_vs_ref64', 'abs_err_max'])),
        'DoR_significant_digits_min': min(
            col(nz, ['DoR_fp32_vs_ref64', 'significant_digits_agreeing_l2'])),
        'DD_significant_digits_min': min(
            [d for d in col(nz, ['DD_fp32_vs_ref64', 'significant_digits_agreeing_l2'])
             if isinstance(d, float)] or [float('inf')]),
        'DD_exact_conditions': sum(
            1 for r in nz
            if r['DD_fp32_vs_ref64']['significant_digits_agreeing_l2'] == 'exact'),
        'downstream_DoR_vs_DD_max_logit_diff': max(
            col(nz, ['downstream', 'DoR_vs_DD', 'max_abs_logit_diff'])),
        'downstream_DD_vs_ref64_max_logit_diff': max(
            col(nz, ['downstream', 'DD_vs_ref64rounded', 'max_abs_logit_diff'])),
        'downstream_DoR_vs_ref64_max_logit_diff': max(
            col(nz, ['downstream', 'DoR_vs_ref64rounded', 'max_abs_logit_diff'])),
        'downstream_DD_bitwise_equal_to_ref64_count': sum(
            1 for r in nz if r['downstream']['DD_vs_ref64rounded']['bitwise_equal']),
        'downstream_DoR_bitwise_equal_to_ref64_count': sum(
            1 for r in nz if r['downstream']['DoR_vs_ref64rounded']['bitwise_equal']),
        'noop_DoR_all_exactly_zero': all(r['noop']['DoR_delta_exactly_zero'] for r in noop),
        'noop_DD_all_exactly_zero': all(r['noop']['DD_delta_exactly_zero'] for r in noop),
        'delta_confined_DoR_all': all(r['DoR_delta_confined'] for r in results),
        'delta_confined_DD_all': all(r['DD_delta_confined'] for r in results),
        'elapsed_s': round(time.time() - t0, 1),
    }

    report = {
        'purpose': 'Compare difference-of-reconstructions (DoR) against the direct '
                   'decoder-delta (DD) implementation of the error-preserving SAE edit, '
                   'against a float64 reference, on the completed validation inputs.',
        'reference_definition': 'float64 direct decoder-delta (DD64). Justified by the '
                                'float64 DoR-vs-DD cross-check, reported per condition.',
        'relative_error_definitions': {
            'rel_err_l2': '||approx - ref||_2 / ||ref||_2 (primary; undefined only if '
                          'the reference delta is identically zero, reported as null)',
            'rel_err_max_vs_peak': 'max|approx - ref| / max|ref|',
            'elementwise': f'|approx_i - ref_i| / |ref_i| computed only where '
                           f'|ref_i| >= {TAU_REL} * max|ref|; coordinates below that '
                           f'threshold are excluded and counted',
            'significant_digits_agreeing_l2': '-log10(rel_err_l2); derived from the '
                                              'RELATIVE error, never from absolute error',
        },
        'edit_strengths': {'amplify': AMPLIFY_FACTOR, 'ablate': 0.0, 'no_op': 1.0,
                           'note': 'unchanged from the completed validation'},
        'aggregate': agg,
        'conditions': results,
    }
    with open(os.path.join(OUT, 'edit_implementation_comparison.json'), 'w') as f:
        json.dump(report, f, indent=1)

    print('\n' + '=' * 76)
    for k, v in agg.items():
        print(f'  {k:46s} {v}')
    print('=' * 76)
    print(f'-> intervention_validation/out/edit_implementation_comparison.json')


if __name__ == '__main__':
    main()
