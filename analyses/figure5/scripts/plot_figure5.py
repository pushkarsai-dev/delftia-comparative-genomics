"""Render the archived domain assignments; never rerun InterProScan."""
import argparse,json,os,tempfile
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'figure5-mpl-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import font_manager
from figure5_data import load_data,e_label,ROOT

COLORS={'Gene3D':'#c4daeb','SUPERFAMILY':'#d6e5ef','Pfam':'#98bfdb','PANTHER':'#b0cbe0'}
WIDTH_MM,HEIGHT_MM=170,160

def build_figure(proteins):
    # Select only actual TTF files: avoid embedding CFF data as TrueType.
    font=font_manager.findfont(font_manager.FontProperties(family='DejaVu Sans'),fallback_to_default=False)
    if Path(font).suffix.lower()!='.ttf':raise ValueError('A TrueType font is required')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
    fig=plt.figure(figsize=(WIDTH_MM/25.4,HEIGHT_MM/25.4))
    for i,p in enumerate(proteins):
        bottom=.565 if i==0 else .085
        ax=fig.add_axes([.205,bottom,.565,.33])
        ax.set_xlim(.5,329.5);ax.set_ylim(-.65,5.1)
        fig.text(.035,bottom+.37,f'({"ab"[i]})',fontsize=11,weight='bold')
        fig.text(.205,bottom+.37,f'{p["id"]}  ·  {p["length"]} aa',fontsize=10,weight='bold')
        fig.text(.205,bottom+.338,'LOG-family signatures' if i==0 else 'Fe(II)/2-oxoglutarate dioxygenase signatures',fontsize=8.5)
        ax.add_patch(Rectangle((.5,4.26),p['length'],.12,facecolor='#d9d9d9',edgecolor='none'))
        if 'motif' in p:
            m=p['motif'];ax.plot([m['start']-.5,m['end']+.5],[4.32,4.32],lw=2.5,color='#b22222')
            ax.text((m['start']+m['end'])/2,4.55,m['sequence'],ha='center',va='bottom',fontsize=8,color='#9c1c1c')
        for k,h in enumerate(p['hits']):
            y=3-k
            if len(h['fragments'])>1:
                ax.plot([h['start']-.5,h['end']+.5],[y+.28,y+.28],color='#52697b',lw=.7,zorder=1)
            for a,b in h['fragments']:
                ax.add_patch(Rectangle((a-.5,y),b-a+1,.56,facecolor=COLORS[h['db']],edgecolor='none',zorder=2))
            ax.text(-.035,(y+.28+.65)/5.75,h['db'],transform=ax.transAxes,ha='right',va='center',fontsize=8)
            largest=max(h['fragments'],key=lambda f:f[1]-f[0])
            ax.text(sum(largest)/2,y+.28,h['accession'],ha='center',va='center',fontsize=7.5,color='#10243a',zorder=3)
            spans=', '.join(f'{a}–{b}' for a,b in h['fragments'])
            ax.text(1.035,(y+.46+.65)/5.75,spans,transform=ax.transAxes,va='center',fontsize=7.5)
            ax.text(1.035,(y+.10+.65)/5.75,e_label(h['evalue']),transform=ax.transAxes,va='center',fontsize=7.5)
        ticks=[1]+list(range(50,p['length'],50))+[p['length']]
        ax.set_xticks(ticks);ax.tick_params(axis='x',labelsize=8,length=2,pad=3)
        ax.set_yticks([])
        for spine in ('left','right','top'):ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_bounds(.5,p['length']+.5)
        ax.set_xlabel('Residue position',fontsize=8.5,labelpad=5)
    return fig

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir',type=Path,default=ROOT/'data')
    ap.add_argument('--output-dir',required=True,type=Path)
    args=ap.parse_args()
    proteins=load_data(args.data_dir)
    dest=args.output_dir.expanduser().resolve()
    protected=[ROOT/'data',ROOT/'reference',ROOT/'scripts',ROOT/'tests',ROOT/'docs']
    if dest.exists() or any(dest==p or p in dest.parents for p in protected):
        ap.error('Choose a fresh output directory outside protected source folders')
    fig=build_figure(proteins)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();bounds=fig.bbox
    for text in fig.findobj(plt.Text):
        if not text.get_visible() or not text.get_text():continue
        b=text.get_window_extent(renderer)
        if b.x0<0 or b.y0<0 or b.x1>bounds.x1 or b.y1>bounds.y1:raise ValueError('Text outside page: '+text.get_text())
    dest.mkdir(parents=True,exist_ok=False)
    fig.savefig(dest/'figure5.pdf',metadata={'Creator':'Figure 5 reproducibility workflow','CreationDate':None,'ModDate':None})
    fig.savefig(dest/'figure5.png',dpi=600)
    (dest/'plotted_data.json').write_text(json.dumps(proteins,indent=2)+'\n')
    plt.close(fig)
    print('Wrote figure5.pdf, figure5.png, and plotted_data.json; 8 signatures, 9 matched segments.')

if __name__=='__main__':main()
