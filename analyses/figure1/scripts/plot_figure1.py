"""Draw Figure 1 from the archived tree, metadata, and directional ANI matrix.

Inputs resolve relative to this repository, independently of the current working
directory. No phylogenetic inference or ANI calculation is performed. The
distance-preserving 40-tip tree keeps unreported support values missing; see
docs/provenance.md for the source and transformation history.

Run with --output-dir naming a new or empty directory. Importing this module
builds the figure in memory without writing figure or mapping files. artifacts()
exposes the actual plotted artists for independent validation.
"""
import argparse, csv, json, re, sys, hashlib, numpy as np
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.transforms import blended_transform_factory

MM = 1 / 25.4
W_MM, H_MM = 170.0, 220.0

ROOT = Path(__file__).resolve().parents[1]
NWK = ROOT / "data/tree_40tips_corrected.newick"
NWK_SHA256_EXPECTED = "ace5e0da9aa2ebcd4689a95d81d4f3f4f2bb4b77a950603f6c5b5c27e456498a"
TIPS = ROOT / "data/tip_mapping.tsv"
METADATA = ROOT / "data/panel_metadata.tsv"
ANI = ROOT / "data/ani40_matrix.tsv"

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.fontset": "custom", "mathtext.rm": "sans", "mathtext.it": "sans:italic",
    "axes.linewidth": 0.6,
    "figure.constrained_layout.use": True,
    "figure.constrained_layout.h_pad": 0.02,
    "figure.constrained_layout.w_pad": 0.02,
})

# ----------------------------------------------------------------- sources ---
def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()

# ------------------------------------------------------------------ newick ---
class Node:
    __slots__ = ("name", "length", "children", "x", "y", "support", "labelled")
    def __init__(self):
        self.name = None; self.length = 0.0; self.children = []
        self.support = None; self.labelled = False

def parse_newick(s):
    s = s.strip().rstrip(";"); pos = 0
    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                n.children.append(node())
                if s[pos] == ",": pos += 1
                else: break
            assert s[pos] == ")"; pos += 1
        m = re.compile(r"[^,():]*").match(s, pos); lab = m.group(0); pos = m.end()
        n.labelled = lab != ""          # absence of a label means support was NOT reported
        if s[pos:pos+1] == ":":
            pos += 1
            m = re.compile(r"[-\d.eE+]+").match(s, pos)
            n.length = float(m.group(0)); pos = m.end()
        n.name = lab or None
        return n
    return node()

def layout(root):
    leaves = []
    def walk(n, depth):
        n.x = depth + n.length
        if not n.children:
            n.y = len(leaves); leaves.append(n); return n.y, n.y
        ys = [walk(c, n.x) for c in n.children]
        n.y = (min(a for a, _ in ys) + max(b for _, b in ys)) / 2
        return min(a for a, _ in ys), max(b for _, b in ys)
    root.length = 0.0
    walk(root, 0.0)
    return leaves

# -------------------------------------------------------------------- data ---
# Pin every compact input, not just the corrected tree. No files are written.
_checksums = json.loads((ROOT / "data/input_checksums.json").read_text())
for _name, _expected in _checksums.items():
    _path = (ROOT / _name).resolve()
    if not _path.is_relative_to(ROOT / "data"):
        raise SystemExit(f"Invalid input path in manifest: {_name}")
    if sha(_path) != _expected:
        raise SystemExit(f"Input checksum mismatch: {_name}")

STUDY_PID = {"BAA554": "BAA554_hybrid", "ATCC13751_own": "ATCC13751_hybrid",
             "ATCC15668_own": "ATCC15668_hybrid", "ALG19": "ALG19_illumina"}
STUDY_KEYS = set(STUDY_PID.values())

def panel():
    with METADATA.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = [r for r in rows if r.get("Tier") in ("T1", "T2", "T3")]
    assert len(rows) == 40, len(rows)
    out = {}
    for r in rows:
        acc = str(r["Assembly accession"] or ""); pid = str(r["Panel ID"])
        m = re.match(r"(GCA_\d+\.\d+)", acc)
        if pid in STUDY_PID:          key = STUDY_PID[pid]
        elif m:                       key = m.group(1)
        elif "zenodo" in acc.lower(): key = "ALG19_illumina"
        else:                         key = pid
        g = str(r["GTDB call (fill)"] or "") + " " + str(r["Species (NCBI)"] or "")
        strain = str(r["Strain"])
        if strain.strip().lower() in ("na", "n/a", "none", ""):
            strain = str(r["Panel ID"])
        if key in out:
            raise ValueError(f"Duplicate panel key: {key}")
        out[key] = dict(strain=strain, panel_id=pid, accession=acc,
                        species="D. acidovorans" if "acidovorans" in g else "D. tsuruhatensis",
                        tier=r["Tier"], type_material=str(r.get("Type material") or "").strip(),
                        study=key in STUDY_KEYS)
    return out

def tip_to_key():
    out = {}
    for r in csv.DictReader(open(TIPS), delimiter="\t"):
        if r["in_panel"] != "yes":
            continue
        lab = r["label"]
        m = re.match(r"(GCA_\d+\.\d+)", lab)
        out[r["tip_id"]] = m.group(1) if m else lab.split(".fna")[0].split(" (User")[0]
    return out

P, T2K = panel(), tip_to_key()
with open(NWK) as _fh:
    _nwk_text = _fh.read()
_nwk_sha = sha(NWK)
if _nwk_sha != NWK_SHA256_EXPECTED:
    raise SystemExit(f"tree source hash mismatch\n  expected {NWK_SHA256_EXPECTED}\n  got      {_nwk_sha}")
root = parse_newick(_nwk_text)
leaves = layout(root)
assert len(leaves) == 40, len(leaves)
keys = [T2K[l.name] for l in leaves]
N = len(leaves)

# Short display labels. All 40 assemblies stay separate; nothing is merged.
SHORT = {"ATCC 15668 (= 2167 = DSM 39 = NBRC 14950)": "ATCC 15668",
         "ATCC BAA-554 (= T7, type strain)": "BAA-554"}
def short(k):
    return SHORT.get(P[k]["strain"], P[k]["strain"])

# Documented culture-collection deposit groups. Each is evidenced by a catalogue
# record, not by ANI. These are the groups DOCUMENTED here; the panel has not
# been exhaustively audited for every possible deposit equivalence.
DEPOSIT_GROUP = {
    # dagger: D. tsuruhatensis T7 type strain, three deposits.
    #   BCCM/LMG catalogue for LMG 29104: Type strain Yes; Depositor No. NBRC
    #   16741; other collection numbers NBRC 16741, DSM 17581, ATCC BAA-554.
    "BAA554_hybrid": "\u2020", "GCA_001571325.1": "\u2020", "GCA_900111975.1": "\u2020",
    # double dagger: D. acidovorans type strain.
    #   Panel Strain field: "ATCC 15668 (= 2167 = DSM 39 = NBRC 14950)".
    "ATCC15668_hybrid": "\u2021", "GCA_059084005.1": "\u2021",
    # section: D. lacustris type strain.
    #   DSMZ catalogue for DSM 21246 gives strain 332 and LMG 24775.
    "GCA_040267515.1": "\u00a7", "GCA_900107225.1": "\u00a7",
}
DEPOSIT_TAXON = {"\u2020": "D. tsuruhatensis", "\u2021": "D. acidovorans",
                 "\u00a7": "D. lacustris"}

# ---------------------------------------------------------------- ANI data ---
rows_ = list(csv.reader(open(ANI), delimiter="\t"))
names = rows_[0][1:]
idx = {n: i for i, n in enumerate(names)}
raw = {r[0]: r[1:] for r in rows_[1:]}

def aniname(k):
    """Accession first, so a strain label can never shadow another genome."""
    for cand in (k, SHORT.get(P[k]["strain"], P[k]["strain"]), P[k]["strain"]):
        if cand in idx:
            return cand
    raise KeyError(f"{k} not in the ANI matrix")

order = [aniname(k) for k in keys]
assert len(set(order)) == N and set(order) == set(names)
M = np.array([[float(raw[q][idx[r]]) for r in order] for q in order])

SPCOL = {"D. tsuruhatensis": "#1F4E9C", "D. acidovorans": "#A0522D"}
SPMARK = {"D. tsuruhatensis": "s", "D. acidovorans": "D"}
norm = Normalize(vmin=95.0, vmax=100.0)
cmap = plt.get_cmap("cividis")          # perceptually uniform, CVD-safe

# ------------------------------------------------------------------ figure ---
fig = plt.figure(figsize=(W_MM * MM, H_MM * MM))
gs = fig.add_gridspec(2, 3, width_ratios=[0.86, 0.72, 1.55],
                      height_ratios=[1.0, 0.085], wspace=0.02, hspace=0.045)
axT = fig.add_subplot(gs[0, 0])
axL = fig.add_subplot(gs[0, 1])
axA = fig.add_subplot(gs[0, 2])
axS = fig.add_subplot(gs[1, 0], sharex=axT)      # scale bar rides axT's x-scale
fig.add_subplot(gs[1, 1]).axis("off")
fig.add_subplot(gs[1, 2]).axis("off")

# ---- (a) tree
SUPPORT_AUDIT = []      # (clade tuple, class, value or None)
n_filled = n_open = n_low = n_unreported = n_root = 0
def tips_of(n):
    return tuple(sorted(_tips(n)))
def _tips(n):
    if not n.children: return [n.name]
    out = []
    for c in n.children: out.extend(_tips(c))
    return out
def draw(n):
    global n_filled, n_open, n_low, n_unreported, n_root
    if not n.children:
        return
    ys = [c.y for c in n.children]
    axT.plot([n.x, n.x], [min(ys), max(ys)], color="#333333", lw=.8, solid_capstyle="butt")
    for c in n.children:
        axT.plot([n.x, c.x], [c.y, c.y], color="#333333", lw=.8, solid_capstyle="butt")
        draw(c)
    if n is root:
        SUPPORT_AUDIT.append((tips_of(n), "root_no_value", None)); n_root += 1
        return
    if not n.labelled:
        # No support value in the inference. NOT low support, and NOT 1.0.
        SUPPORT_AUDIT.append((tips_of(n), "unreported", None)); n_unreported += 1
        return
    sup = float(n.name)
    SUPPORT_AUDIT.append((tips_of(n), "measured", sup))
    if sup >= 0.95:
        axT.plot(n.x, n.y, "o", ms=3.0, mfc="#111111", mec="none", zorder=5); n_filled += 1
    elif sup >= 0.70:
        axT.plot(n.x, n.y, "o", ms=3.0, mfc="white", mec="#111111", mew=.6, zorder=5); n_open += 1
    else:
        n_low += 1
draw(root)

xmax = max(l.x for l in leaves)
for l in leaves:
    axT.plot([l.x, xmax * 1.02], [l.y, l.y], color="#D6D6D6", lw=.4, zorder=0)
axT.set_xlim(-xmax * 0.03, xmax * 1.03)
axT.set_ylim(N - 0.5, -0.5)
axT.axis("off")
axT.text(-0.01, 1.004, "(a)", transform=axT.transAxes, fontsize=10,
         fontweight="bold", va="bottom", ha="left")

# ---- scale bar: data x on axT's shared scale, axes y on its own small axis
SB = round(xmax * 0.25, 3)
axS.set_ylim(0, 1)
axS.axis("off")
bt = blended_transform_factory(axS.transData, axS.transAxes)
x_lo = 0.0
axS.plot([x_lo, x_lo + SB], [0.78, 0.78], transform=bt, color="#333333", lw=1.0,
         solid_capstyle="butt", clip_on=False)
for xe in (x_lo, x_lo + SB):
    axS.plot([xe, xe], [0.66, 0.90], transform=bt, color="#333333", lw=1.0, clip_on=False)
axS.text(x_lo + SB / 2, 0.42, f"{SB:g} substitutions per site", transform=bt,
         ha="center", va="top", fontsize=8)

# ---- shared numbered labels
axL.set_xlim(0, 1); axL.set_ylim(N - 0.5, -0.5); axL.axis("off")
for i, (l, k) in enumerate(zip(leaves, keys), start=1):
    d = P[k]
    axL.plot(0.055, l.y, marker=SPMARK[d["species"]], ms=3.2,
             mfc=SPCOL[d["species"]], mec=SPCOL[d["species"]], mew=.5, clip_on=False)
    txt = f"{i}  {short(k)}"
    if k in DEPOSIT_GROUP:
        txt += " " + DEPOSIT_GROUP[k]
    axL.text(0.14, l.y, txt, ha="left", va="center", fontsize=8,
             color=SPCOL[d["species"]], fontweight="bold" if d["study"] else "normal")
    if d["study"]:
        axL.plot(0.985, l.y, marker="*", ms=5.0, mfc="#B22222", mec="none", clip_on=False)

# ---- (b) directional ANI
im = axA.imshow(M, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest",
                extent=[0, N, N - 0.5, -0.5])
for v in range(0, N + 1, 5):                      # guides every 5, not per cell
    axA.plot([v, v], [-0.5, N - 0.5], color="white", lw=.45, alpha=.75, zorder=3)
    axA.plot([0, N], [v - 0.5, v - 0.5], color="white", lw=.45, alpha=.75, zorder=3)
axA.set_xlim(0, N); axA.set_ylim(N - 0.5, -0.5)
axA.set_yticks([])
# Sparse column labels. A tick marks EVERY column so each one is locatable;
# only the numbers in LABELLED_COLUMNS are printed, which removes the crowding
# of forty labels across a 78 mm axis. Row numbers beside the tree give the
# mapping for every column.
LABELLED_COLUMNS = [1] + list(range(5, N + 1, 5))
axA.set_xticks(np.arange(N) + .5)
axA.set_xticklabels([])
axA.tick_params(axis="x", length=1.6, width=.5, color="#999999", pad=1.0)
for i in LABELLED_COLUMNS:
    k = keys[i - 1]
    axA.plot([i - .5, i - .5], [0, -0.012], transform=blended_transform_factory(
        axA.transData, axA.transAxes), color="#555555", lw=.7, clip_on=False, zorder=6)
    axA.text(i - .5, -0.016, str(i),
             transform=blended_transform_factory(axA.transData, axA.transAxes),
             ha="center", va="top", fontsize=8, clip_on=False,
             color=SPCOL[P[k]["species"]],
             fontweight="bold" if P[k]["study"] else "normal")
for s_ in axA.spines.values():
    s_.set_linewidth(.6); s_.set_color("#888888")
axA.set_xlabel("reference genome (row number)", fontsize=9, labelpad=16)
axA.text(-0.028, 0.5, "query genome (row number)", transform=axA.transAxes,
         rotation=90, ha="center", va="center", fontsize=9)
axA.text(-0.045, 1.004, "(b)", transform=axA.transAxes, fontsize=10,
         fontweight="bold", va="bottom", ha="left")

cb = fig.colorbar(im, ax=axA, fraction=.030, pad=.015)
cb.set_label("directional ANI (%)", fontsize=9, labelpad=3)
cb.set_ticks([95, 96, 97, 98, 99, 100])
cb.ax.tick_params(labelsize=8, length=2.0, width=.6)
cb.outline.set_linewidth(.6)

# ---- key
def _rows_for(sym):
    return [i for i, k in enumerate(keys, start=1) if DEPOSIT_GROUP.get(k) == sym]

def _dep_legend(sym):
    """Row numbers and the deposit count are read from the drawn order."""
    r = _rows_for(sym)
    assert len(r) >= 2, f"{sym} marks {len(r)} rows"
    # Oxford comma: "38, 39, and 40" for three or more, "30 and 31" for two.
    joined = (f"{r[0]} and {r[1]}" if len(r) == 2
              else ", ".join(str(x) for x in r[:-1]) + f", and {r[-1]}")
    return f"{sym} rows {joined}: one strain, {len(r)} deposits"

h = [plt.Line2D([], [], marker=SPMARK["D. tsuruhatensis"], ms=4.4, ls="none",
                mfc=SPCOL["D. tsuruhatensis"], mec=SPCOL["D. tsuruhatensis"]),
     plt.Line2D([], [], marker=SPMARK["D. acidovorans"], ms=4.4, ls="none",
                mfc=SPCOL["D. acidovorans"], mec=SPCOL["D. acidovorans"]),
     plt.Line2D([], [], marker="*", ms=7, mfc="#B22222", mec="none", ls="none"),
     plt.Line2D([], [], marker="o", ms=4.0, mfc="#111111", mec="none", ls="none"),
     plt.Line2D([], [], marker="o", ms=4.0, mfc="white", mec="#111111", mew=.7, ls="none"),
     plt.Line2D([], [], ls="none", marker="none"),
     plt.Line2D([], [], ls="none", marker="none"),
     plt.Line2D([], [], ls="none", marker="none"),
     plt.Line2D([], [], ls="none", marker="none")]
labels = ["D. tsuruhatensis", "D. acidovorans",
          "sequenced or reassembled in this study",
          "local support \u2265 0.95",
          "local support 0.70 \u2264 s < 0.95",
          "unmarked: support < 0.70, or not reported",
          _dep_legend("\u2020"), _dep_legend("\u2021"), _dep_legend("\u00a7")]
leg = fig.legend(h, labels, loc="lower left", bbox_to_anchor=(0.235, 0.002),
                 frameon=False, fontsize=8, handlelength=1.1, labelspacing=.40,
                 ncols=2, columnspacing=1.1, borderaxespad=0.0)
for t in list(leg.get_texts())[:2]:
    t.set_style("italic")

def artifacts():
    """Live objects for validation.

    The validator inspects THESE, not a re-derivation of the inputs, so a
    transpose or a mis-set extent in the plotting code is detectable.
    """
    return {
        "fig": fig, "root": root, "ax_tree": axT, "ax_labels": axL, "ax_ani": axA, "ax_scale": axS,
        "image": im,                      # the heatmap artist itself
        "keys": list(keys),               # row order actually drawn
        "order": list(order),             # ANI matrix names, row order
        "leaves": list(leaves),
        "tree_sha256": _nwk_sha,
        "tree_path": NWK,
        "scale_bar": SB, "xmax": xmax,
        "support_audit": list(SUPPORT_AUDIT),
        "deposit_group": dict(DEPOSIT_GROUP),
        "labelled_columns": list(LABELLED_COLUMNS),
        "legend_texts": [t.get_text() for t in leg.get_texts()],
        "species_colour": dict(SPCOL), "species_marker": dict(SPMARK),
        "norm": norm, "cmap_name": cmap.name if hasattr(cmap, "name") else str(cmap),
        "panel": P,
    }


def _write_outputs(output_dir):
    output_dir = Path(output_dir).expanduser().resolve()
    for protected in (ROOT / "data", ROOT / "reference", ROOT / "scripts", ROOT / "docs", ROOT / "tests"):
        if output_dir == protected or protected in output_dir.parents:
            raise ValueError(f"Output directory is inside protected repository content: {output_dir}")
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        raise ValueError(f"Output directory must be new or empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "figure1.pdf")
    fig.savefig(output_dir / "figure1.png", dpi=600)
    with (output_dir / "support_audit.csv").open("w", newline="", encoding="utf-8") as _fh:
        _w = csv.writer(_fh)
        _w.writerow(["n_tips", "support_class", "support_value", "tips"])
        for cl, kind, val in sorted(SUPPORT_AUDIT, key=lambda r: (len(r[0]), r[0])):
            _w.writerow([len(cl), kind, "" if val is None else f"{val:g}", " ".join(cl)])
    mapping = {
            "tree_source": NWK.relative_to(ROOT).as_posix(),
            "tree_source_sha256": _nwk_sha,
            "scale_bar_substitutions_per_site": SB,
            "tree_xmax": xmax,
            "labelled_columns": LABELLED_COLUMNS,
            "deposit_groups": {sym: {"rows": _rows_for(sym),
                                     "taxon": DEPOSIT_TAXON[sym],
                                     "n_deposits": len(_rows_for(sym))}
                               for sym in sorted(set(DEPOSIT_GROUP.values()))},
            "support_counts": {"measured_ge_0.95": n_filled,
                               "measured_0.70_to_0.95": n_open,
                               "measured_lt_0.70": n_low,
                               "unreported": n_unreported,
                               "root_no_value": n_root},
            "row_numbers": {str(i): {"label": short(k), "key": k,
                                     "panel_id": P[k]["panel_id"],
                                     "accession": P[k]["accession"],
                                     "ani_matrix_name": order[i-1],
                                     "tip_id": leaves[i-1].name,
                                     "species": P[k]["species"],
                                     "study": P[k]["study"],
                                     "deposit_group": DEPOSIT_GROUP.get(k)}
                            for i, k in enumerate(keys, start=1)},
        }
    with (output_dir / "row_mapping.json").open("w", encoding="utf-8") as _out:
        json.dump(mapping, _out, indent=1, ensure_ascii=False)
        _out.write("\n")
    with (output_dir / "row_mapping.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["row_number", "display_label", "species", "study", "panel_id",
                         "assembly_accession", "ani_matrix_name", "kbase_tip_id",
                         "deposit_symbol", "column_number_printed"])
        for i in range(1, N + 1):
            r = mapping["row_numbers"][str(i)]
            writer.writerow([i, r["label"], r["species"], "yes" if r["study"] else "no",
                             r["panel_id"], r["accession"], r["ani_matrix_name"], r["tip_id"],
                             r["deposit_group"] or "", "yes" if i in LABELLED_COLUMNS else "no"])
    print(f"markers: filled {n_filled}  open {n_open}  unmarked-low {n_low}  "
          f"unmarked-unreported {n_unreported}  root {n_root}")
    print(f"figure {fig.get_size_inches()[0]*25.4:.2f} x {fig.get_size_inches()[1]*25.4:.2f} mm")
    print("wrote figure1.pdf/.png, row_mapping.json/.csv, and support_audit.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build/figure1",
                        help="New or empty output directory (default: module/build/figure1)")
    args = parser.parse_args()
    try:
        _write_outputs(args.output_dir)
    except (OSError, ValueError) as exc:
        parser.exit(2, f"Error: {exc}\n")
