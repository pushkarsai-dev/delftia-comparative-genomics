"""Render Figure 2 from independently recomputed, checksum-verified inputs."""
import argparse
import json
import platform
import sys
from pathlib import Path
import matplotlib
from figstyle import W, panel_label, save, font
import matplotlib.pyplot as plt
from figure2_data import ROOT, check_destination, compute, write_tables

PAN, CORE = "#1F5A92", "#B25A1E"
CLS = [("cloud",     "#F0C89A", "cloud, 1–5 genomes"),
       ("shell",     "#D98E45", "shell, 6–37"),
       ("soft_core", "#5B8CC0", "soft core, 38–39"),
       ("core",      "#1F5A92", "core, all 40")]


def create_figure(result):
    nod = result['clades']
    acc = result['means']
    gain, loss = result['gain'], result['core_loss']
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(W, 105/25.4),
                                   gridspec_kw=dict(width_ratios=[1, 1.08], wspace=0.12))

    # ---------------------------------------------------------------- panel (a) ---
    n   = [int(r["n_genomes"]) for r in acc]
    pan = [r["pan_mean"] for r in acc]
    cm  = [r["core_mean"] for r in acc]

    axA.scatter([int(r["n_genomes"]) for r in nod], [int(r["total_families"]) for r in nod],
                s=12, fc="none", ec=PAN, lw=.65, zorder=2)
    axA.scatter([int(r["n_genomes"]) for r in nod], [int(r["perfect_core"]) for r in nod],
                s=12, fc="none", ec=CORE, lw=.65, zorder=2)
    # Individual-clade points display the full observed spread; no range band.
    axA.plot(n, pan, "-o", color=PAN,  ms=3.2, lw=1.0, zorder=3, label="pangenome")
    axA.plot(n, cm,  "-o", color=CORE, ms=3.2, lw=1.0, zorder=3, label="core genome")

    axA.annotate(f"{pan[-1]:,}", (40, pan[-1]), textcoords="offset points",
                 xytext=(-2, 6), ha="right", fontsize=8, color=PAN, fontweight="bold")
    axA.annotate(f"{cm[-1]:,}", (40, cm[-1]), textcoords="offset points",
                 xytext=(-2, -12), ha="right", fontsize=8, color=CORE, fontweight="bold")
    axA.annotate(f"39 → 40 genomes:\n+{gain:,} families,\n−{loss:,} from core", (40, 12000),
                 textcoords="offset points", xytext=(-4, 0), ha="right", va="center",
                 fontsize=8, color="#444444", linespacing=1.3)

    axA.set_xlabel("genomes in clade", fontsize=9)
    axA.set_ylabel("orthologous families", fontsize=9)
    axA.set_xlim(0, 43); axA.set_ylim(0, 23600)
    axA.set_xticks([2, 10, 20, 30, 40])
    axA.set_yticks(range(0, 24000, 5000))
    axA.set_yticklabels([f"{v:,}" for v in range(0, 24000, 5000)])
    for s in ("top", "right"): axA.spines[s].set_visible(False)
    axA.legend(loc="upper left", frameon=False, fontsize=8, handlelength=1.5,
               borderpad=0, labelspacing=.3)
    panel_label(axA, "a", x=-0.20, y=1.00)

    # ---------------------------------------------------------------- panel (b) ---
    cnt = {n: result["occupancy"][str(n)] for n in range(1, 41)}
    cls = {}
    for g in cnt: cls[g] = 'core' if g == 40 else 'soft_core' if g >= 38 else 'shell' if g >= 6 else 'cloud'
    xs = list(range(1, 41))
    cols = [dict((c, col) for c, col, _ in CLS)[cls[g]] for g in xs]
    axB.bar(xs, [cnt[g] for g in xs], color=cols, width=.82, lw=0)
    axB.set_yscale("log")
    axB.set_ylim(8, 40000)
    axB.set_xlim(0.2, 40.8)
    axB.set_xticks([1, 10, 20, 30, 40])
    axB.set_xlabel("genomes carrying the family", fontsize=9)
    axB.set_ylabel("families (log scale)", fontsize=9)
    for s in ("top", "right"): axB.spines[s].set_visible(False)

    axB.annotate(f"{cnt[1]:,}", (1, cnt[1]), textcoords="offset points", xytext=(1, 3),
                 ha="left", fontsize=8, color="#8A5A20", fontweight="bold")
    axB.annotate(f"{cnt[40]:,}", (40, cnt[40]), textcoords="offset points", xytext=(-1, 3),
                 ha="right", fontsize=8, color=PAN, fontweight="bold")

    h = [plt.Rectangle((0, 0), 1, 1, fc=col) for _, col, _ in CLS]
    tot = result['categories']
    lab = [f"{t}   {tot[k]:,}" for k, _, t in CLS]
    axB.legend(h[::-1], lab[::-1], loc="upper center", frameon=False, fontsize=8,
               handlelength=1.0, handleheight=.9, borderpad=0, labelspacing=.32,
               bbox_to_anchor=(0.52, 1.02))
    panel_label(axB, "b", x=-0.155, y=1.00)

    return fig

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'build/figure2')
    args = parser.parse_args()
    try:
        out = check_destination(args.output_dir)
        result = compute()
        fig = create_figure(result)
        out.mkdir(parents=True, exist_ok=False)
        write_tables(result, out)
        save(fig, 'figure2', out)
        (out/'execution.json').write_text(json.dumps(dict(python=platform.python_version(),
            matplotlib=matplotlib.__version__, font=font, module='figure2',
            input_sha256=result['input_sha256']), indent=2)+'\n')
        plt.close(fig)
        print(f'Created {out / "figure2.pdf"}')
    except (ValueError, OSError, KeyError, TypeError, StopIteration) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    return 0

if __name__ == '__main__':
    sys.exit(main())
