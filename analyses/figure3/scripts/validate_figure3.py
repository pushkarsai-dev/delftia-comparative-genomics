"""Fail on cell, label, export, font, or reproducibility mismatches."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import matplotlib.pyplot as plt
import numpy as np
import pdfplumber
from PIL import Image
import pypdfium2

from figure3_data import FAMILIES, MODULE, STUDY, EXCEPTION, load_data
from plot_figure3 import build_figure, COLORS, BOUNDS

def validate_artists(fig,ax,genomes,values):
    expected=np.asarray(values)
    if not np.array_equal(ax.images[0].get_array(),expected):
        raise ValueError('plotted matrix differs')
    cells={(round(t.get_position()[1]-.5),round(t.get_position()[0]-.5)):t
           for t in ax.texts if t.get_position()[1]<40}
    if len(cells)!=440:
        raise ValueError('expected 440 cell labels')
    for i in range(40):
        for j in range(11):
            if cells[i,j].get_text()!=str(values[i][j]):
                raise ValueError('cell label differs')
    expected_labels=[STUDY.get(g,g)+(' *' if g in STUDY else ' †' if g==EXCEPTION else '') for g in genomes]
    if [t.get_text() for t in ax.get_yticklabels()]!=expected_labels:
        raise ValueError('genome labels differ')
    if [t.get_text() for t in ax.get_xticklabels()]!=[f'{f}  {label}' for f,label,_ in FAMILIES]:
        raise ValueError('family labels differ')
    for text,g in zip(ax.get_yticklabels(),genomes):
        if (text.get_fontweight()=='bold') != (g in STUDY or g==EXCEPTION):
            raise ValueError('genome highlighting differs')
    if list(ax.images[0].cmap.colors)!=COLORS or not np.array_equal(ax.images[0].norm.boundaries,BOUNDS):
        raise ValueError('color scale differs')
    if not np.allclose(fig.get_size_inches()*25.4,[170,220]):
        raise ValueError('page size differs')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for text in [*ax.texts,*ax.get_xticklabels(),*ax.get_yticklabels()]:
        bb=text.get_window_extent(renderer)
        if bb.x0<0 or bb.y0<0 or bb.x1>fig.bbox.width or bb.y1>fig.bbox.height:
            raise ValueError(f'clipped figure text: {text.get_text()}')
    if any(t.get_fontsize()<8 for t in [*ax.texts,*ax.get_xticklabels(),*ax.get_yticklabels()]):
        raise ValueError('essential text below 8 pt')

def raster_pdf(path):
    doc=pypdfium2.PdfDocument(str(path))
    return np.array(doc[0].render(scale=2).to_pil())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--pdffonts',default='pdffonts')
    args=p.parse_args()
    genomes,values=load_data()
    fig,ax=build_figure(genomes,values)
    validate_artists(fig,ax,genomes,values)
    out=args.output_dir
    with (out/'displayed_counts.tsv').open() as handle:
        rows=list(csv.reader(handle,delimiter='\t'))
    assert rows==[['genome']+[f[0] for f in FAMILIES]]+[[g]+list(map(str,v)) for g,v in zip(genomes,values)]
    # Independently count raw hit records instead of trusting the rerun summary.
    observed=Counter()
    with (MODULE/'data/independent_hits.tsv').open() as handle:
        for hit in csv.DictReader(handle,delimiter='\t'):
            if hit['mode']=='per_genome' and float(hit['evalue'])<1e-3:
                observed[hit['target'].split('|')[0],hit['family']]+=1
    for i,g in enumerate(genomes):
        for j,(f,_,_) in enumerate(FAMILIES):
            assert observed[g,f]==values[i][j]
    with pdfplumber.open(out/'figure3.pdf') as doc:
        assert len(doc.pages)==1
        page=doc.pages[0]
        assert abs(page.width*25.4/72-170)<.001 and abs(page.height*25.4/72-220)<.001
        assert not [c for c in page.chars if c['x0']<-.01 or c['x1']>page.width+.01 or c['top']<-.01 or c['bottom']>page.height+.01]
        # Extract the exported glyphs inside the known heatmap bounds.
        left,bottom,width,height=.275,.095,.707,.70
        glyphs=[c for c in page.chars if c['upright'] and c['text'].isdigit()
                and left*page.width<c['x0']<(left+width)*page.width
                and (1-bottom-height)*page.height<c['top']<(1-bottom)*page.height]
        assert len(glyphs)==440, len(glyphs)
        extracted={}
        for c in glyphs:
            j=int(((c['x0']+c['x1'])/2/page.width-left)/width*11)
            i=int(((c['top']+c['bottom'])/2/page.height-(1-bottom-height))/height*40)
            assert (i,j) not in extracted
            extracted[i,j]=int(c['text'])
        assert all(extracted[i,j]==values[i][j] for i in range(40) for j in range(11))
    fonts=subprocess.run([args.pdffonts,str(out/'figure3.pdf')],capture_output=True,text=True,check=True)
    assert not fonts.stderr.strip(),fonts.stderr
    assert 'TrueType' in fonts.stdout and 'Type 0C' not in fonts.stdout
    with tempfile.TemporaryDirectory(prefix='figure3-export-check-') as tmp:
        tmp=Path(tmp)
        fig.savefig(tmp/'check.pdf')
        fig.savefig(tmp/'check.png',dpi=600)
        assert np.array_equal(np.asarray(Image.open(out/'figure3.png')),np.asarray(Image.open(tmp/'check.png')))
        assert np.array_equal(raster_pdf(out/'figure3.pdf'),raster_pdf(tmp/'check.pdf'))
    plt.close(fig)
    print('PASS: 440 matrix cells, labels, raw-hit counts, and exported PDF numbers; 170 x 220 mm; no clipped text; embedded TrueType; fresh PDF/PNG renders identical.')

if __name__=='__main__':
    main()
