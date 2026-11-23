#!/usr/bin/env python3
"""Figure: per-layer P3 AUROC. Reads ONLY layer_sweep/v1/p3_layer_curves.tsv."""
import csv, sys, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

V = sys.argv[1]
r = list(csv.DictReader(open(os.path.join(V, 'p3_layer_curves.tsv')), delimiter='\t'))
L = sorted({int(x['layer']) for x in r})
g = lambda fs: [float([z for z in r if int(z['layer']) == l and z['feature_set'] == fs][0]['auroc_P3']) for l in L]
comp, esm, both = g('comp'), g('esm'), g('comp_plus_esm')
SEL, REF = 26, 8

THEMES = {
 'light': dict(surface='#fcfcfb', ink='#0b0b0b', ink2='#52514e', grid='#e3e2de',
               s1='#2a78d6', s2='#eb6834', s3='#1baf7a'),
 'dark':  dict(surface='#1a1a19', ink='#f0efec', ink2='#c3c2b7', grid='#333330',
               s1='#3987e5', s2='#d95926', s3='#199e70'),
}
for mode, t in THEMES.items():
    fig, ax = plt.subplots(figsize=(9.2, 5.0), dpi=200)
    fig.patch.set_facecolor(t['surface']); ax.set_facecolor(t['surface'])
    ax.grid(True, color=t['grid'], lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(t['grid']); ax.spines[sp].set_linewidth(1.0)

    ax.plot(L, comp, color=t['s3'], lw=2.0, ls=(0, (5, 3)), zorder=3,
            label='composition only (41 features, layer-independent)')
    ax.plot(L, esm, color=t['s1'], lw=2.0, marker='o', ms=4.0, zorder=5,
            markeredgecolor=t['surface'], markeredgewidth=1.0, label='ESM only (2,560)')
    ax.plot(L, both, color=t['s2'], lw=2.0, marker='s', ms=3.6, zorder=4,
            markeredgecolor=t['surface'], markeredgewidth=1.0,
            label='composition + ESM (2,601)')

    for x, c in ((REF, t['ink2']), (SEL, t['ink2'])):
        ax.axvline(x, color=c, lw=1.0, ls=':', zorder=1, alpha=0.7)
    i8, i26 = L.index(REF), L.index(SEL)
    ax.annotate(f'layer 8\n{both[i8]:.4f}', (REF, both[i8]), textcoords='offset points',
                xytext=(0, -34), ha='center', fontsize=8.5, color=t['ink2'])
    ax.annotate(f'layer 26 — selected\n{both[i26]:.4f}', (SEL, both[i26]),
                textcoords='offset points', xytext=(0, 14), ha='center',
                fontsize=8.5, color=t['ink'], fontweight='bold')

    # direct labels (relief for the contrast WARN on the aqua slot)
    ax.annotate('composition only', (L[-1], comp[-1]), textcoords='offset points',
                xytext=(8, 0), va='center', fontsize=8.5, color=t['s3'])
    ax.annotate('ESM only', (L[-1], esm[-1]), textcoords='offset points',
                xytext=(8, 5), va='center', fontsize=8.5, color=t['s1'])
    ax.annotate('comp + ESM', (L[-1], both[-1]), textcoords='offset points',
                xytext=(8, -7), va='center', fontsize=8.5, color=t['s2'])

    ax.set_xlabel('ESM-2-650M encoder layer  (hidden_states[i] = output of block i)',
                  fontsize=9.5, color=t['ink2'])
    ax.set_ylabel('AUROC on P3 (n = 2,548)', fontsize=9.5, color=t['ink2'])
    ax.set_title('Localisation decodability by layer — probes fitted on P1+P2, scored on P3',
                 fontsize=11.5, color=t['ink'], pad=12, loc='left')
    ax.set_xticks([1, 4, 8, 12, 16, 20, 24, 26, 29, 33])
    ax.tick_params(colors=t['ink2'], labelsize=8.5)
    ax.set_xlim(0, 39); ax.set_ylim(0.940, 0.9925)
    leg = ax.legend(loc='lower right', frameon=False, fontsize=8.5, ncol=1)
    for txt in leg.get_texts(): txt.set_color(t['ink2'])
    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(V, f'fig_p3_layer_curve_{mode}.{ext}'),
                    facecolor=t['surface'])
    plt.close(fig)
    print(f'wrote fig_p3_layer_curve_{mode}.png/.svg')
