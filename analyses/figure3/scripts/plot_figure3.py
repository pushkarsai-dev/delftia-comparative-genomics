"""Render Figure 3 without modifying archived data or accepted references."""
import argparse
import csv
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.font_manager import FontProperties, findfont
import numpy as np

from figure3_data import FAMILIES, STUDY, EXCEPTION, load_data, reserve_output

COLORS = ['#F0F0F0','#D6E6F4','#ABCDE8','#78ABD5','#4682BA','#1F5A92']
BOUNDS = [0,1,2,3,5,7,100]

def build_figure(genomes, values):
    font = 'DejaVu Sans'
    for candidate in ['Arial','Liberation Sans','DejaVu Sans']:
        try:
            path = findfont(FontProperties(family=candidate),fallback_to_default=False)
            if Path(path).suffix.lower() == '.ttf':
                font = candidate
                break
        except ValueError:
            pass
    plt.rcParams.update({'font.family':font,'font.size':8,'pdf.fonttype':42,
                         'ps.fonttype':42,'figure.constrained_layout.use':False})
    fig = plt.figure(figsize=(170/25.4,220/25.4))
    ax = fig.add_axes([0.275,0.095,0.707,0.70])
    matrix = np.asarray(values,dtype=int)
    im = ax.imshow(matrix,cmap=ListedColormap(COLORS),norm=BoundaryNorm(BOUNDS,6),
                   interpolation='nearest',aspect='auto',extent=[0,11,40,0])
    for i in range(40):
        for j in range(11):
            v = matrix[i,j]
            ax.text(j+.5,i+.5,str(v),ha='center',va='center',fontsize=8,
                    color='white' if v>=5 else '#555555' if v==0 else '#171717')
    ax.set_xticks(np.arange(11)+.5)
    ax.set_xticklabels([f'{f}  {label}' for f,label,_ in FAMILIES],rotation=90,
                       ha='left',va='bottom',rotation_mode='anchor',fontsize=8)
    ax.tick_params(axis='x',labeltop=True,labelbottom=False,pad=4,length=0)
    ax.set_yticks(np.arange(40)+.5)
    ax.set_yticklabels([STUDY.get(g,g)+(' *' if g in STUDY else ' †' if g==EXCEPTION else '') for g in genomes],fontsize=8)
    ax.tick_params(axis='y',length=0,pad=4)
    for text,g in zip(ax.get_yticklabels(),genomes):
        if g in STUDY:
            text.set_fontweight('bold')
        if g==EXCEPTION:
            text.set_color('#A52620');text.set_fontweight('bold')
    for x in range(12):
        ax.axvline(x,color='white',lw=.55)
    for y in range(41):
        ax.axhline(y,color='white',lw=.5)
    for x in [2,4,6,9]:
        ax.axvline(x,color='#555555',lw=.8)
    ax.axhline(4,color='#555555',lw=.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Compact grouping key shares the heatmap's column coordinates.
    groups = [(0,2,'Targets'),(2,4,'UbiD/UbiX\nrelated'),(4,6,'Other\nDCs'),(6,9,'Trp\nmetabolism'),(9,11,'Controls')]
    for start,end,label in groups:
        ax.plot([start+.12,end-.12],[40.35,40.35],color='#555555',lw=.7,clip_on=False)
        ax.text((start+end)/2,40.55,label,ha='center',va='top',fontsize=8,clip_on=False)
    ax.set_xlim(0,11)
    ax.set_ylim(40,0)
    cax = fig.add_axes([0.34,0.025,0.49,0.012])
    cb = fig.colorbar(im,cax=cax,orientation='horizontal',spacing='uniform',ticks=[.5,1.5,2.5,4,6,53.5])
    cb.ax.set_xticklabels(['0','1','2','3–4','5–6','≥7'],fontsize=8)
    cb.ax.tick_params(length=0,pad=2)
    cb.outline.set_linewidth(.5)
    fig.text(.585,.043,'Detected proteins per genome',ha='center',fontsize=8)
    return fig,ax

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,required=True)
    args = p.parse_args()
    try:
        genomes,values = load_data()
        out = reserve_output(args.output_dir)
    except (ValueError,OSError) as e:
        p.error(str(e))
    fig,ax = build_figure(genomes,values)
    fig.savefig(out/'figure3.pdf',metadata={'Title':'Figure 3. Detected protein-family counts','Author':''})
    fig.savefig(out/'figure3.png',dpi=600)
    with (out/'displayed_counts.tsv').open('w') as handle:
        writer=csv.writer(handle,delimiter='\t'); writer.writerow(['genome']+[f[0] for f in FAMILIES])
        writer.writerows([[g]+v for g,v in zip(genomes,values)])
    (out/'execution.json').write_text(json.dumps(dict(python=platform.python_version(),matplotlib=matplotlib.__version__,numpy=np.__version__,
        width_mm=170,height_mm=220,orientation='genomes rows, families columns',cells=440),indent=2)+'\n')
    print(f'Wrote Figure 3, 440 cells, and execution record to {out}')

if __name__=='__main__':
    main()
