#!/usr/bin/env python3
"""Read-only data, artist, export, and caption validation for Figure 1.

Run after ``python scripts/plot_figure1.py``. Paths are repository-relative;
the current working directory does not matter. The plotter is imported to
inspect its actual artists, but importing it does not export files. Tests
include negative controls for transposition, one altered ANI cell, and a
wrongly filled support marker. No validation failure regenerates an output.
"""

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

# Importing the plotting module must not create repository bytecode artifacts.
sys.dont_write_bytecode = True

import numpy as np
import pdfplumber
import pypdfium2
from PIL import Image

REPO = Path(__file__).resolve().parents[1]
ANI_TSV = REPO / "data/ani40_matrix.tsv"
CAPTION = REPO / "docs/figure1.md"
TREE_SHA_EXPECTED = "ace5e0da9aa2ebcd4689a95d81d4f3f4f2bb4b77a950603f6c5b5c27e456498a"
APPROVED_FONTS = ("Arial", "Helvetica", "LiberationSans", "DejaVuSans")
FAIL = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAIL.append(name)
    return bool(ok)


def expected_matrix(order):
    """Independent header-keyed lookup, not the plotter's positional index."""
    with ANI_TSV.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        first = reader.fieldnames[0]
        rows = list(reader)
        by_query = {row[first]: row for row in rows}
    if len(by_query) != len(rows):
        raise ValueError("ANI input contains duplicate query identifiers")
    return np.array([[float(by_query[q][r]) for r in order] for q in order])


def compare_plotted(plotted, expected, tol=0.0):
    if plotted.shape != expected.shape:
        return expected.size, float("inf")
    diff = np.abs(np.asarray(plotted, dtype=float) - expected)
    finite = np.isfinite(diff)
    return int(((diff > tol) | ~finite).sum()), float(diff.max())


def walk(node):
    yield node
    for child in node.children:
        yield from walk(child)


def expected_marker_classes(root):
    result = {}
    for node in walk(root):
        if not node.children or node is root:
            continue
        pos = (round(node.x, 9), round(node.y, 9))
        if not node.labelled:
            result[pos] = "none"
        else:
            value = float(node.name)
            result[pos] = "filled" if value >= 0.95 else "open" if value >= 0.70 else "none"
    return result


def plotted_marker_classes(axis):
    result = {}
    for line in axis.get_lines():
        if line.get_marker() != "o" or line.get_linestyle() not in ("None", "none", ""):
            continue
        xs, ys = line.get_xdata(), line.get_ydata()
        if len(xs) != 1:
            continue
        face = str(line.get_markerfacecolor()).lower()
        kind = "open" if face in ("white", "w", "#ffffff") else "filled"
        pos = (round(float(xs[0]), 9), round(float(ys[0]), 9))
        if pos in result:
            raise ValueError(f"Duplicate support marker at {pos}")
        result[pos] = kind
    return result


def compare_markers(expected, plotted):
    bad = [(pos, want, plotted.get(pos, "none")) for pos, want in expected.items()
           if plotted.get(pos, "none") != want]
    bad.extend((pos, "none", plotted[pos]) for pos in plotted if pos not in expected)
    return bad


def validate_mapping(artifacts, mapping, output_dir):
    numbers = mapping["row_numbers"]
    expected_rows = {str(i) for i in range(1, 41)}
    check("mapping has exactly rows 1 through 40", set(numbers) == expected_rows)
    if set(numbers) != expected_rows:
        return
    rows = [numbers[str(i)] for i in range(1, 41)]
    check("mapping row order equals the drawn key order", [r["key"] for r in rows] == artifacts["keys"])
    check("mapping ANI names equal the plotted column order",
          [r["ani_matrix_name"] for r in rows] == artifacts["order"])
    check("mapping tip order equals the plotted leaf order",
          [r["tip_id"] for r in rows] == [leaf.name for leaf in artifacts["leaves"]])
    check("40 unique tips, keys, panel IDs, and matrix names",
          all(len({r[field] for r in rows}) == 40 for field in ("tip_id", "key", "panel_id", "ani_matrix_name")))
    check("mapping tree digest matches the pinned tree", mapping["tree_source_sha256"] == TREE_SHA_EXPECTED)
    check("mapping scale bar equals the drawn scale bar",
          mapping["scale_bar_substitutions_per_site"] == artifacts["scale_bar"])
    check("mapping column-label scheme equals the drawn scheme",
          mapping["labelled_columns"] == artifacts["labelled_columns"])
    check("mapping deposit symbols equal the plotter's assignments",
          all(r.get("deposit_group") == artifacts["deposit_group"].get(r["key"]) for r in rows))
    check("mapping metadata agrees with the panel actually used for plotting",
          all(all(r[field] == artifacts["panel"][r["key"]][field]
                  for field in ("panel_id", "accession", "species", "study")) for r in rows))
    audit = artifacts["support_audit"]
    values = [value for _, kind, value in audit if kind == "measured"]
    expected_counts = {"measured_ge_0.95": sum(value >= .95 for value in values),
                       "measured_0.70_to_0.95": sum(.70 <= value < .95 for value in values),
                       "measured_lt_0.70": sum(value < .70 for value in values),
                       "unreported": sum(kind == "unreported" for _, kind, _ in audit),
                       "root_no_value": sum(kind == "root_no_value" for _, kind, _ in audit)}
    check("saved support counts agree with the plotted clade audit",
          mapping.get("support_counts") == expected_counts
          and all(type(value) is int for value in mapping.get("support_counts", {}).values()),
          str(expected_counts))
    with (output_dir / "row_mapping.csv").open(newline="", encoding="utf-8") as fh:
        csv_rows = list(csv.DictReader(fh))
    # The human-readable CSV uses descriptive column names for JSON fields.
    check("mapping CSV has 40 rows", len(csv_rows) == 40)
    fields = {"ani_matrix_name": "ani_matrix_name", "kbase_tip_id": "tip_id",
              "panel_id": "panel_id", "species": "species", "display_label": "label",
              "assembly_accession": "accession"}
    consistent = len(csv_rows) == 40 and all(
        all(c.get(csv_field) == r[json_field] for csv_field, json_field in fields.items())
        and c.get("deposit_symbol") == (r.get("deposit_group") or "")
        and c.get("study") == ("yes" if r["study"] else "no")
        for c, r in zip(csv_rows, rows))
    check("mapping CSV identifiers and order agree with JSON", consistent)
    check("mapping CSV row numbers are 1 through 40",
          [r.get("row_number") for r in csv_rows] == [str(i) for i in range(1, 41)])
    check("mapping CSV printed-column flags agree with actual scheme",
          all(r.get("column_number_printed") == ("yes" if i in artifacts["labelled_columns"] else "no")
              for i, r in enumerate(csv_rows, start=1)))


def decoded_png(data):
    """Compare image content, not PNG compression or metadata."""
    with Image.open(io.BytesIO(data)) as image:
        return np.asarray(image.convert("RGBA")).copy()


def rendered_pdf(data, dpi=200):
    """Render every PDF page entirely in memory, without writing temporary files."""
    images = []
    with pypdfium2.PdfDocument(data) as document:
        for page in document:
            try:
                bitmap = page.render(scale=dpi / 72)
                try:
                    images.append(np.asarray(bitmap.to_pil().convert("RGBA")).copy())
                finally:
                    bitmap.close()
            finally:
                page.close()
    return images


def pixel_difference(actual, expected):
    """Return differing-pixel count, or None if geometry differs."""
    if actual.shape != expected.shape:
        return None
    return int(np.any(actual != expected, axis=-1).sum())


def export_correspondence(artifacts, pdf_path, png_path):
    """Bind saved exports to the artists checked against scientific inputs.

    Generate fresh exports in memory in the plotter's exact PDF-then-PNG order.
    Full-image comparisons require generation and validation in the same
    dependency/font environment; metadata timestamps do not affect the checks.
    """
    print("\n6A. SAVED EXPORTS vs FRESH IN-MEMORY EXPORTS")
    fresh_pdf, fresh_png = io.BytesIO(), io.BytesIO()
    artifacts["fig"].savefig(fresh_pdf, format="pdf")
    artifacts["fig"].savefig(fresh_png, format="png", dpi=600)
    actual = decoded_png(png_path.read_bytes())
    expected = decoded_png(fresh_png.getvalue())
    difference = pixel_difference(actual, expected)
    check("all saved PNG pixels equal the current figure export", difference == 0,
          f"{difference} differing pixels" if difference is not None else "image dimensions differ")
    del actual, expected
    actual_pages = rendered_pdf(pdf_path.read_bytes())
    expected_pages = rendered_pdf(fresh_pdf.getvalue())
    if check("saved PDF page count equals the current figure export", len(actual_pages) == len(expected_pages)):
        for index, (actual, expected) in enumerate(zip(actual_pages, expected_pages), start=1):
            difference = pixel_difference(actual, expected)
            check(f"all rendered PDF pixels equal current export, page {index} (200 dpi)", difference == 0,
                  f"{difference} differing pixels" if difference is not None else "page dimensions differ")


def validate_exports(artifacts, pdf_path, png_path, pdffonts):
    """Inspect the saved PDF/PNG independently of in-memory geometry."""
    bar_length = None
    with pdfplumber.open(pdf_path) as document:
        if not check("PDF contains exactly one page", len(document.pages) == 1):
            return
        page = document.pages[0]
        width, height = page.width, page.height
        print("\n6. SPARSE COLUMN LABELS")
        wanted = [1] + list(range(5, 41, 5))
        check("labelled columns are 1, 5, 10 ... 40", artifacts["labelled_columns"] == wanted,
              str(artifacts["labelled_columns"]))
        band = [c for c in page.chars if c["upright"] and c["text"].isdigit()
                and c["x0"] > 235 and 520 < c["top"] < 560]
        runs, current = [], []
        for char in sorted(band, key=lambda c: c["x0"]):
            if current and char["x0"] - current[-1]["x1"] > 1.2:
                runs.append(current)
                current = []
            current.append(char)
        if current:
            runs.append(current)
        got = ["".join(c["text"] for c in run) for run in runs]
        check("exactly nine specified numbers are printed under the heatmap", got == list(map(str, wanted)), str(got))
        gaps = [runs[i + 1][0]["x0"] - runs[i][-1]["x1"] for i in range(len(runs) - 1)]
        check("no collision between printed column numbers", bool(gaps) and min(gaps) > 0,
              f"min gap {min(gaps):.2f} pt = {min(gaps) / 72 * 25.4:.2f} mm" if gaps else "no measurable gaps")
        short_vertical = [line for line in page.lines if abs(line["x0"] - line["x1"]) < .3
                          and abs(line["y1"] - line["y0"]) < 4 and line["x0"] > 235]
        bands = Counter(round(line["top"], 1) for line in short_vertical)
        top, count = bands.most_common(1)[0] if bands else (None, 0)
        check("a tick is drawn for all 40 columns", count == 40, f"{count} ticks at top={top}")

        print("\n7. PAGE GEOMETRY AND TEXT")
        check("page is 170.00 x 220.00 mm",
              abs(width / 72 * 25.4 - 170) < .05 and abs(height / 72 * 25.4 - 220) < .05,
              f"{width / 72 * 25.4:.2f} x {height / 72 * 25.4:.2f} mm")
        horizontal = [line for line in page.lines if abs(line["y0"] - line["y1"]) < .3]
        vertical = [line for line in page.lines if abs(line["x0"] - line["x1"]) < .3]
        bars = sorted([line for line in horizontal if line["x1"] - line["x0"] > 5], key=lambda line: line["y0"])
        if check("a PDF scale-bar candidate is present", bool(bars)):
            bar = bars[0]
            bar_length = bar["x1"] - bar["x0"]
            th = [line for line in horizontal if line["x1"] < 175 and line["y0"] > bar["y0"] + 8]
            tv = [line for line in vertical if line["x0"] < 175 and line["y0"] > bar["y0"] + 8]
            if check("independent tree-scale anchors are present", bool(th) and bool(tv)):
                scale = (max(line["x1"] for line in th) - min(line["x0"] for line in tv)) / (artifacts["xmax"] * 1.02)
                expected = artifacts["scale_bar"] * scale
                check("scale bar length matches the tree scale", abs(bar_length - expected) < 1e-6,
                      f"{bar_length:.4f} pt, expected {expected:.4f} pt, error {100 * (bar_length - expected) / expected:+.4f}%")
        objects = list(page.chars) + list(page.lines) + [rect for rect in page.rects
                   if not (rect["width"] > width * .98 and rect["height"] > height * .98)]
        if check("PDF has measurable text and graphic objects", bool(objects)):
            xs = [obj[key] for obj in objects for key in ("x0", "x1")]
            ys = [obj[key] for obj in objects for key in ("top", "bottom")]
            check("nothing clipped", min(xs) > 0 and max(xs) < width and min(ys) > 0 and max(ys) < height,
                  f"margins mm L{min(xs) / 72 * 25.4:.2f} R{(width - max(xs)) / 72 * 25.4:.2f} "
                  f"T{min(ys) / 72 * 25.4:.2f} B{(height - max(ys)) / 72 * 25.4:.2f}")
        upright = sorted({round(char["size"], 1) for char in page.chars if char["upright"]})
        check("upright text is 7 pt or larger", bool(upright) and min(upright) >= 7.0, str(upright))
        groups = defaultdict(list)
        for char in page.chars:
            if not char["upright"]:
                groups[round(char["x0"], 0)].append(char)
        effective = [round(max(c["x1"] for c in chars) - min(c["x0"] for c in chars), 2)
                     for chars in groups.values() if len(chars) >= 4]
        check("rotated labels meet the 8 pt minimum (x-extent at 90 degrees)",
              bool(effective) and min(effective) >= 8.0, f"measured sizes {sorted(effective)} pt")

    print("\n8. FONTS")
    result = subprocess.run([pdffonts, str(pdf_path)], capture_output=True, text=True, check=False)
    check("pdffonts exited 0", result.returncode == 0, f"rc={result.returncode}")
    check("no font-type mismatch warning", "mismatch" not in result.stderr.lower(), result.stderr.strip() or "clean stderr")
    lines = [line for line in result.stdout.splitlines()[2:] if line.strip()]
    check("at least one font is listed", bool(lines), f"{len(lines)} fonts")
    pattern = re.compile(r"^(\S+)\s+(.+?)\s+(\S+)\s+(yes|no)\s+(yes|no)\s+(yes|no)\s+\d+\s+\d+\s*$")
    parsed = [pattern.match(line) for line in lines]
    check("font table is parseable", bool(parsed) and all(parsed))
    approved = bool(parsed) and all(match and "TrueType" in match.group(2)
                                   and match.group(4) == "yes"
                                   and any(face in match.group(1).replace(" ", "") for face in APPROVED_FONTS)
                                   for match in parsed)
    check("every font is an approved face, explicitly embedded, and TrueType", approved,
          "; ".join(line.split()[0] for line in lines))

    print("\n9. PNG")
    with Image.open(png_path) as image:
        pixels_w, pixels_h = image.size
        array = np.asarray(image.convert("L"), dtype=float)
    check("PNG is 600 dpi at the same physical size",
          abs(pixels_w / 600 * 25.4 - 170) < .1 and abs(pixels_h / 600 * 25.4 - 220) < .1,
          f"{pixels_w}x{pixels_h}px = {pixels_w / 600 * 25.4:.2f} x {pixels_h / 600 * 25.4:.2f} mm")
    pixels_per_point = pixels_w / width
    center_length = None
    for y in range(pixels_h - 1, int(pixels_h * .80), -1):
        dark = np.where(array[y, :int(pixels_w * .40)] < 120)[0]
        if len(dark) > 50 and dark.max() - dark.min() == len(dark) - 1:
            center_length = dark.max() - dark.min() + 1 - pixels_per_point
            break
    if check("PNG scale bar is independently detected", center_length is not None):
        if check("PDF scale-bar measurement is available for PNG comparison", bar_length is not None):
            diff_pixels = abs(center_length - bar_length * pixels_per_point)
            check("PNG scale bar agrees with PDF within one pixel at 600 dpi", diff_pixels <= 1.0,
                  f"{diff_pixels:.2f} px = {diff_pixels / pixels_per_point / 72 * 25.4:.4f} mm")


def validate_caption(mapping):
    print("\n10. CANONICAL CAPTION AND ALT TEXT (docs/figure1.md)")
    blocks = re.findall(r"```(?:text)?\n(.*?)\n```", CAPTION.read_text(encoding="utf-8"), re.S)
    if not check("caption, alt text, and Methods blocks are present", len(blocks) >= 3, f"{len(blocks)} fenced blocks"):
        return
    caption, alt = blocks[:2]
    groups_match = True
    for symbol, group in sorted(mapping["deposit_groups"].items()):
        rows = group["rows"]
        if len(rows) == 2:
            phrase = f"{symbol} rows {rows[0]} and {rows[1]} are two deposits of one strain"
        elif len(rows) == 3:
            phrase = f"{symbol} rows {rows[0]}, {rows[1]}, and {rows[2]} are three deposits of one strain"
        else:
            groups_match = False
            continue
        groups_match = groups_match and phrase in caption
    check("caption deposit groups match rendered mapping exactly", groups_match)
    check("caption states missing and low support separately",
          "Internal nodes without a circle either have support < 0.70 or lack a reported support value; the root also has no support value." in caption)
    check("caption describes the exported sparse column scheme",
          "a tick marks every column and selected column numbers are labeled" in caption)
    check("caption scale-bar value matches mapping",
          f"{mapping['scale_bar_substitutions_per_site']:g} substitutions per site" in caption)
    check("caption limits deposit claim rather than asserting an exhaustive audit", "exhaustive" in caption)
    check("alt text is at most 255 characters", len(alt) <= 255, f"{len(alt)} characters")
    check("alt text says 'sequenced or reassembled'", "sequenced or reassembled" in alt)
    check("caption paragraphs have no internal line breaks", all("\n" not in para for para in caption.split("\n\n")))
    stale = [text for text in ("rows 26 and 27", "rows 16 and 18", "38, 39 and 40") if text in caption]
    check("no superseded row list survives in caption", not stale, str(stale))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=REPO / "build/figure1", help="Existing figure-output directory")
    parser.add_argument("--pdffonts", default="pdffonts", help="Poppler pdffonts executable or absolute path")
    args = parser.parse_args(argv)
    output_dir = args.output_dir.resolve()
    required = [ANI_TSV, CAPTION, output_dir / "figure1.pdf", output_dir / "figure1.png",
                output_dir / "row_mapping.json", output_dir / "row_mapping.csv", output_dir / "support_audit.csv"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        print("Cannot validate: required files are missing:\n  " + "\n  ".join(missing), file=sys.stderr)
        return 2
    pdffonts = shutil.which(args.pdffonts)
    if pdffonts is None:
        print("Cannot validate: Poppler pdffonts is unavailable; install Poppler or use --pdffonts PATH.", file=sys.stderr)
        return 2
    FAIL.clear()
    spec = importlib.util.spec_from_file_location("plot_figure1", Path(__file__).with_name("plot_figure1.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    artifacts = module.artifacts()
    print("=" * 74 + "\n1. TREE SOURCE")
    check("tree digest equals pinned expected digest", artifacts["tree_sha256"] == TREE_SHA_EXPECTED)
    digest = hashlib.sha256(Path(artifacts["tree_path"]).read_bytes()).hexdigest()
    check("digest independently recomputed from disk matches", digest == TREE_SHA_EXPECTED)
    check("40 tips drawn", len(artifacts["leaves"]) == 40, str(len(artifacts["leaves"])))

    print("\n2. HEATMAP ARTIST vs INDEPENDENT EXPECTED MATRIX")
    image = artifacts["image"]
    plotted = np.asarray(image.get_array(), dtype=float)
    expected = expected_matrix(artifacts["order"])
    bad, worst = compare_plotted(plotted, expected)
    check("plotted array equals expected matrix in every cell", bad == 0,
          f"{plotted.size} cells, {bad} mismatched, worst |diff| {worst:g}")
    transposed_bad, _ = compare_plotted(plotted, expected.T)
    check("NEGATIVE CONTROL: rejects transposed matrix", transposed_bad > 0, f"{transposed_bad} cells mismatch")
    corrupt = expected.copy()
    corrupt[7, 19] += .01
    altered_bad, _ = compare_plotted(plotted, corrupt)
    check("NEGATIVE CONTROL: rejects one altered cell", altered_bad == 1)
    extent = tuple(float(value) for value in image.get_extent())
    check("extent maps row i to y=i and column j to x=j+0.5", extent == (0., 40., 39.5, -.5), str(extent))
    check("y axis puts row 1 at top", artifacts["ax_ani"].get_ylim() == (39.5, -.5))
    check("normalization is 95 to 100", (image.norm.vmin, image.norm.vmax) == (95., 100.))
    check("no cell clipped by normalization", plotted.min() >= 95 and plotted.max() <= 100,
          f"{plotted.min():.4f}-{plotted.max():.4f}")
    check("matrix not symmetrized", np.abs(plotted - plotted.T).max() > .2,
          f"max |M-M.T| {np.abs(plotted - plotted.T).max():.4f}")
    check("diagonal not forced to 100", np.diag(plotted).min() < 100, f"min {np.diag(plotted).min():.4f}")

    print("\n3. ROW / COLUMN MAPPING")
    mapping = json.loads((output_dir / "row_mapping.json").read_text(encoding="utf-8"))
    validate_mapping(artifacts, mapping, output_dir)
    print("\n4. SUPPORT MARKERS, READ FROM PLOTTED ARTISTS")
    expected_classes = expected_marker_classes(artifacts["root"])
    plotted_classes = plotted_marker_classes(artifacts["ax_tree"])
    bad_markers = compare_markers(expected_classes, plotted_classes)
    check("every plotted marker matches tree support label", not bad_markers,
          f"{len(plotted_classes)} markers, {len(bad_markers)} disagreements")
    counts = Counter(expected_classes.values())
    check("classes are 3 filled, 17 open, 15 unmarked",
          (counts["filled"], counts["open"], counts["none"]) == (3, 17, 15), str(dict(counts)))
    wrong = dict(expected_classes)
    wrong[next(pos for pos, kind in wrong.items() if kind == "none")] = "filled"
    check("NEGATIVE CONTROL: wrongly filled marker rejected", len(compare_markers(wrong, plotted_classes)) == 1)
    audit = artifacts["support_audit"]
    unreported = [row for row in audit if row[1] == "unreported"]
    check("10 clades have unreported support and no numeric value", len(unreported) == 10 and all(row[2] is None for row in unreported))
    expected_csv = [(str(len(clade)), kind, "" if value is None else f"{value:g}", " ".join(clade))
                    for clade, kind, value in sorted(audit, key=lambda row: (len(row[0]), row[0]))]
    with (output_dir / "support_audit.csv").open(newline="", encoding="utf-8") as fh:
        saved_audit = list(csv.DictReader(fh))
    actual_csv = [(row["n_tips"], row["support_class"], row["support_value"], row["tips"]) for row in saved_audit]
    check("saved support audit agrees with plotted clades", actual_csv == expected_csv)

    print("\n5. LEGEND vs RENDERED MAPPING")
    legend = " ".join(artifacts["legend_texts"])
    groups_match = True
    for symbol, group in sorted(mapping["deposit_groups"].items()):
        rows = group["rows"]
        joined = f"{rows[0]} and {rows[1]}" if len(rows) == 2 else ", ".join(map(str, rows[:-1])) + f", and {rows[-1]}"
        groups_match = groups_match and f"{symbol} rows {joined}: one strain, {len(rows)} deposits" in legend
    check("every deposit group states its rendered rows and count", groups_match)
    check("three-deposit group present and says '3 deposits'", any(len(g["rows"]) == 3 for g in mapping["deposit_groups"].values()) and "3 deposits" in legend)
    check("legend explains low and unreported support", any("support < 0.70, or not reported" in text for text in artifacts["legend_texts"]))
    export_correspondence(artifacts, output_dir / "figure1.pdf", output_dir / "figure1.png")
    validate_exports(artifacts, output_dir / "figure1.pdf", output_dir / "figure1.png", pdffonts)
    validate_caption(mapping)
    print("\n" + "=" * 74)
    if FAIL:
        print(f"FAILED: {len(FAIL)} check(s): {FAIL}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, IndexError, TypeError, RuntimeError) as error:
        print(f"Validation could not complete: {type(error).__name__}: {error}", file=sys.stderr)
        sys.exit(2)
