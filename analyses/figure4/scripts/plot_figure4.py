"""Generate Figure 4 from archived inputs without overwriting existing outputs."""
import argparse, csv, json, platform
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap,Normalize
from figure4_data import ROOT,MARK,load_data,arrow_vertices

NAVY='#16305C';ORANGE='#C9772F';RED='#B22222'
GRAD=LinearSegmentedColormap.from_list('gene_order',['#F4CDA0','#C9772F','#4A7FB5','#16305C'])

def build():
    genes,tracks,totals=load_data()
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans','DejaVu Sans'],
                         'font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
    fig=plt.figure(figsize=(170/25.4,220/25.4),layout=None)
    key=[Patch(color=NAVY),Patch(color=ORANGE),Line2D([],[],marker='|',color='#555555',ls='',markersize=7),Patch(facecolor='none',edgecolor=RED)]
    fig.legend(key,['Detected in all 40 genomes','Detected in 25–39 genomes','BGC0000984 BLAST match','Highlighted protein'],
               loc='upper center',bbox_to_anchor=(.53,.997),frameon=False,ncols=2,fontsize=8,columnspacing=1.2,handlelength=1.2)
    fig.text(.025,.922,'(a)',weight='bold',fontsize=10)
    fig.text(.10,.925,'BAA-554 delftibactin-associated region',weight='bold',fontsize=9)
    axa=fig.add_axes([.10,.833,.875,.062])
    artists=[]
    for g in genes:
        # Source positions are 1-based inclusive; display uses genomic interval edges.
        x0=(g['start']-1)/1000;x1=g['end']/1000
        patch=Polygon(arrow_vertices(x0,x1,0,g['strand'],height=.44,head_limit=.65),closed=True,
                      facecolor=NAVY if g['conservation']==40 else ORANGE,
                      edgecolor=RED if g['gene'] in MARK else 'white',linewidth=.8 if g['gene'] in MARK else .25)
        axa.add_patch(patch);artists.append(('a',g,patch))
        if g['kcb']:axa.plot([(x0+x1)/2]*2,[-.32,-.48],color='#555555',lw=.7)
        if g['gene'] in MARK:
            right=g['order']>30
            axa.annotate(MARK[g['gene']],((x0+x1)/2,.23),xytext=(-10 if right else 0,10),textcoords='offset points',
                         ha='right' if right else 'center',fontsize=8,color=RED,arrowprops={'arrowstyle':'-','color':RED,'lw':.65})
    axa.set_xlim(2158.5,2246.5);axa.set_ylim(-.55,.65);axa.set_yticks([])
    axa.set_xticks([2160,2180,2200,2220,2240],['2.16','2.18','2.20','2.22','2.24'])
    axa.set_xlabel('Chromosome position (Mb)',labelpad=3)
    for s in ['top','right','left']:axa.spines[s].set_visible(False)
    axa.tick_params(labelsize=8,length=3)
    fig.text(.025,.766,'(b)',weight='bold',fontsize=10)
    fig.text(.10,.769,f'Gene order in {len(tracks)} selected assemblies',weight='bold',fontsize=9)
    fig.text(.028,.738,'Genome',fontsize=8,weight='bold')
    fig.text(.838,.738,'Hits',fontsize=8,weight='bold',ha='center')
    fig.text(.949,.738,'Span (kb)',fontsize=8,weight='bold',ha='center')
    axb=fig.add_axes([.295,.32,.505,.405])
    labelax=fig.add_axes([.02,.32,.269,.405],sharey=axb)
    valueax=fig.add_axes([.812,.32,.173,.405],sharey=axb)
    for i,t in enumerate(tracks):
        y=-i
        labelax.text(.965,y,t['label'],ha='right',va='center',fontsize=8,weight='bold' if t['genome'] in ['BAA554_hybrid','ATCC13751_hybrid'] else 'normal')
        if t['reversed']:axb.plot(-1.35,y,marker='<',color=RED,ms=4,clip_on=False)
        for g in t['genes']:
            patch=Polygon(arrow_vertices(g['x0']/1000,g['x1']/1000,y,g['display_strand'],height=.47,head_limit=.8),closed=True,
                          facecolor=GRAD((g['ref_gene_order']-1)/44),edgecolor='white',linewidth=.20)
            axb.add_patch(patch);artists.append(('b',dict(g,genome=t['genome']),patch))
        valueax.text(.15,y,f"{len(t['genes'])}/45",ha='center',va='center',fontsize=8)
        valueax.text(.79,y,f"{t['span_bp']/1000:.1f}",ha='center',va='center',fontsize=8)
    axb.set_xlim(0,86);axb.set_ylim(-len(tracks)+.3,.7)
    labelax.set_xlim(0,1);valueax.set_xlim(0,1)
    for ax in [axb,labelax,valueax]:ax.axis('off')
    # A physical scale is explicit; rows are not independently rescaled.
    scale=fig.add_axes([.295,.300,.505,.01],sharex=axb)
    scale.plot([0,10],[.5,.5],color='black',lw=1);scale.plot([0,0],[.15,.85],color='black',lw=.7);scale.plot([10,10],[.15,.85],color='black',lw=.7)
    scale.text(12,.5,'10 kb',va='center',fontsize=8);scale.set_ylim(0,1);scale.axis('off')
    cax=fig.add_axes([.38,.269,.34,.009])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(1,45),cmap=GRAD),cax=cax,orientation='horizontal',ticks=[1,23,45])
    cb.set_label('Gene position in the BAA-554 region',fontsize=8,labelpad=2);cb.ax.tick_params(labelsize=8,length=2,pad=2);cb.outline.set_linewidth(.5)
    fig.text(.025,.215,'(c)',weight='bold',fontsize=10)
    fig.text(.10,.218,'Protein conservation across all 40 genomes',weight='bold',fontsize=9)
    axc=fig.add_axes([.10,.071,.875,.13])
    bars=axc.bar([g['order'] for g in genes],[g['conservation'] for g in genes],width=.72,color=[NAVY if g['conservation']==40 else ORANGE for g in genes])
    for g in genes:
        if g['conservation']<40:axc.text(g['order'],g['conservation']+.7,str(g['conservation']),ha='center',fontsize=8,color='#8A5A20')
        if g['gene'] in MARK:axc.plot(g['order'],42.5,'v',color=RED,ms=4)
    axc.set_xlim(.3,45.7);axc.set_ylim(0,46);axc.set_yticks([0,20,40]);axc.set_xticks([1,10,20,30,40,45])
    axc.set_xlabel('Gene position in BAA-554 chromosome order',labelpad=3)
    axc.set_ylabel('Genomes with\nan ortholog',labelpad=3)
    axc.tick_params(labelsize=8,length=3)
    for s in ['top','right']:axc.spines[s].set_visible(False)
    audit=dict(genes=genes,tracks=tracks,region_start_1based=2159932,region_end_1based=2245340,
               coordinate_convention='1-based inclusive inputs; 0-based interval edges for arrow geometry; span=end-start+1',
               figure_mm=[170,220],scale_bar_kb=10,whole_proteome_totals=totals)
    return fig,audit,artists,bars,axb,scale

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path);a=p.parse_args()
    if a.output_dir.exists():p.error('Output directory already exists; choose a fresh directory')
    fig,audit,*_=build()
    a.output_dir.mkdir(parents=True)
    fig.savefig(a.output_dir/'figure4.pdf',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(a.output_dir/'figure4.png',dpi=600)
    (a.output_dir/'plot_data.json').write_text(json.dumps(audit,indent=2)+'\n')
    (a.output_dir/'execution.json').write_text(json.dumps({'python':platform.python_version(),'matplotlib':matplotlib.__version__},indent=2)+'\n')
    print(f"Wrote Figure 4: 45 reference genes, {len(audit['tracks'])} tracks, and 45 conservation bars.")
if __name__=='__main__':main()
