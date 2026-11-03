"""Validate Figure 2 values, exported geometry, fonts, and complete rendered images."""
import argparse
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import numpy as np
import pdfplumber
import pypdfium2
from PIL import Image
from matplotlib.colors import to_rgba
from figure2_data import ROOT, category, compute, require, rows
from plot_figure2 import create_figure, CLS
import matplotlib.pyplot as plt

def validate_artists(fig, result):
    a, b = fig.axes
    require(len(a.collections) == 2, 'Expected two scatter series and no shading')
    require(len(a.lines) == 2, 'Expected two mean lines')
    for artist, key in zip(a.collections, ('total_families', 'perfect_core')):
        expected = [(r['n_genomes'], r[key]) for r in result['clades']]
        require(np.array_equal(artist.get_offsets(), expected), 'Scatter values differ')
    for line, key in zip(a.lines, ('pan_mean', 'core_mean')):
        require(np.array_equal(line.get_xdata(), [r['n_genomes'] for r in result['means']]), 'Mean x values differ')
        require(np.array_equal(line.get_ydata(), [r[key] for r in result['means']]), 'Mean y values differ')
    require(len(b.patches) == 40, 'Expected 40 histogram bars')
    colors = {key: value for key, value, _ in CLS}
    for n, bar in enumerate(b.patches, 1):
        require(abs(bar.get_x() + bar.get_width()/2 - n) < 1e-10, 'Bar position differs')
        require(bar.get_height() == result['occupancy'][str(n)], 'Bar height differs')
        require(np.allclose(bar.get_facecolor(), to_rgba(colors[category(n)])), 'Bar color differs')
    require(b.get_yscale() == 'log', 'Histogram axis is not logarithmic')
    require(b.get_ylim()[0] < min(result['occupancy'].values()) and b.get_ylim()[1] > max(result['occupancy'].values()), 'Histogram values clipped')
    require(np.allclose(fig.get_size_inches()*25.4, [170,105]), 'Figure dimensions differ')
    legend = [t.get_text() for t in b.get_legend().get_texts()]
    require(legend == [f'{label}   {result["categories"][key]:,}' for key,_,label in CLS][::-1], 'Category legend differs')
    labels = [t.get_text() for ax in fig.axes for t in ax.texts]
    for value in ('21,853', '4,263', '11,544', '39 → 40 genomes:\n+312 families,\n−3 from core'):
        require(value in labels, f'Missing annotation: {value}')

def validate_exported_bars(bars, axis, result):
    require(len(bars) == 40, 'Expected 40 exported histogram bars')
    colors = {k:to_rgba(v)[:3] for k,v,_ in CLS}
    worst = 0
    for n, r in enumerate(bars, 1):
        x = axis['x0'] + (n-.2)/40.6*(axis['x1']-axis['x0'])
        y = axis['bottom'] - (math.log10(result['occupancy'][str(n)])-math.log10(8))/(math.log10(40000)-math.log10(8))*(axis['bottom']-axis['top'])
        error = max(abs((r['x0']+r['x1'])/2-x), abs(r['top']-y))
        require(error < 1e-4, 'Exported bar geometry differs')
        require(np.allclose(r['non_stroking_color'], colors[category(n)], atol=1e-8, rtol=0), 'Exported color differs')
        worst = max(worst,error)
    return worst

def rendered_pdf(source):
    with pypdfium2.PdfDocument(source) as doc:
        require(len(doc) == 1, 'Expected one PDF page')
        return np.asarray(doc[0].render(scale=2).to_pil().convert('RGB')).copy()

def validate(folder, pdffonts):
    folder = Path(folder)
    result = compute()
    require(json.loads((folder/'summary.json').read_text()) == result, 'Saved summary differs from raw-input recomputation')
    expected_tables = {'clades.tsv': result['clades'], 'means.tsv': result['means'],
                       'occupancy.tsv':[dict(n_genomes=n,n_families=result['occupancy'][str(n)],category=category(n)) for n in range(1,41)]}
    for name, values in expected_tables.items():
        expected = [{k:'; '.join(v) if isinstance(v,list) else str(v) for k,v in row.items()} for row in values]
        require(rows(folder/name) == expected, f'Saved table differs: {name}')
    with pdfplumber.open(folder/'figure2.pdf') as doc:
        require(len(doc.pages) == 1, 'Expected one figure page')
        p = doc.pages[0]
        require(np.allclose([p.width*25.4/72,p.height*25.4/72],[170,105],atol=1e-4), 'PDF dimensions differ')
        bars = sorted((r for r in p.rects if r['height'] > 1000), key=lambda r:r['x0'])
        axes = [r for r in p.rects if r['non_stroking_color'] == 1.0 and r['x0'] > p.width/2]
        require(len(axes) == 1, 'Cannot identify histogram plotting area')
        error = validate_exported_bars(bars,axes[0],result)
        require(all(0 <= c['x0'] <= c['x1'] <= p.width and 0 <= c['top'] <= c['bottom'] <= p.height for c in p.chars), 'PDF text crosses page bounds')
    fonts = subprocess.run([pdffonts,str(folder/'figure2.pdf')], capture_output=True,text=True,check=True)
    require(not fonts.stderr.strip(), 'PDF font warning: '+fonts.stderr)
    font_rows = fonts.stdout.splitlines()[2:]
    require(font_rows and all('TrueType' in r and 'yes' in r for r in font_rows), 'Expected embedded TrueType fonts')
    fig = create_figure(result)
    validate_artists(fig,result)
    pdf, png = io.BytesIO(), io.BytesIO()
    fig.savefig(pdf,format='pdf'); fig.savefig(png,format='png',dpi=600)
    plt.close(fig)
    require(np.array_equal(np.asarray(Image.open(folder/'figure2.png')), np.asarray(Image.open(png))), 'Complete PNG differs from fresh rendering')
    require(np.array_equal(rendered_pdf(str(folder/'figure2.pdf')),rendered_pdf(pdf.getvalue())), 'Complete PDF rendering differs from fresh rendering')
    return dict(status='passed', families=result['families'], matrix_cells=40*result['families'],
                clades=len(result['clades']), scatter_points=72, mean_points=30,
                exported_bars=40, max_bar_geometry_error_points=error,
                full_png_equal=True, full_pdf_render_equal=True, font_warnings=0)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'build/figure2')
    parser.add_argument('--pdffonts', default='pdffonts')
    args = parser.parse_args()
    try:
        print(json.dumps(validate(args.output_dir,args.pdffonts),indent=2))
    except (ValueError,OSError,KeyError,TypeError,StopIteration,subprocess.CalledProcessError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    return 0

if __name__ == '__main__':
    sys.exit(main())
