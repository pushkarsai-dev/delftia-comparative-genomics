"""Figure 2 style, independent of other manuscript figures."""
from pathlib import Path
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import font_manager
W = 170/25.4
font = None
for family in ('Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'):
    try:
        candidate = font_manager.findfont(family, fallback_to_default=False)
    except ValueError:
        continue
    if Path(candidate).suffix.lower() == '.ttf':
        font = family
        break
if font is None: raise RuntimeError('A supported TrueType sans-serif font is required')
mpl.rcParams.update({'font.family':font, 'font.size':8, 'axes.labelsize':9,
    'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42,
    'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,
    'figure.constrained_layout.use':True,'figure.constrained_layout.h_pad':.05,
    'figure.constrained_layout.w_pad':.05,'mathtext.fontset':'dejavusans'})
def panel_label(ax,letter,x=-.10,y=1.00):
    ax.text(x,y,f'({letter})',transform=ax.transAxes,fontsize=10,fontweight='bold',ha='left',va='bottom')
def save(fig,name,outdir):
    outdir=Path(outdir)
    targets=[outdir/(name+'.pdf'),outdir/(name+'.png')]
    if any(p.exists() for p in targets): raise FileExistsError('Refusing to replace existing figure files')
    outdir.mkdir(parents=True,exist_ok=True)
    fig.savefig(targets[0]); fig.savefig(targets[1],dpi=600)
