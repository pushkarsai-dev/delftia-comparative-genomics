"""Cross-check source coordinates, plotted geometry, labels, fonts, and exported files."""
import argparse,csv,json,re,subprocess,tempfile,io
from collections import Counter
from pathlib import Path
import numpy as np
import pdfplumber,pypdfium2
from PIL import Image
from Bio import SeqIO
import matplotlib.pyplot as plt
from figure4_data import ROOT,MARK
from plot_figure4 import build

def independent_sources():
    data=ROOT/'data';ref=SeqIO.read(data/'regions/BAA554_hybrid.gbk','genbank')
    offset=int(ref.annotations['structured_comment']['antiSMASH-Data']['Orig. start'])
    features=sorted([f for f in ref.features if f.type=='CDS'],key=lambda f:int(f.location.start))
    expected={}
    for i,f in enumerate(features,1):
        tag=f.qualifiers['locus_tag'][0].replace('ctg1_','BAA554_hybrid_contig_1_')
        expected[tag]=(offset+int(f.location.start),offset+int(f.location.end),f.location.strand,i)
    order=list(csv.DictReader(io.StringIO((data/'region_gene_order_corrected.tsv').read_text()),delimiter='\t'))
    per={g:[r for r in order if r['genome']==g] for g in {r['genome'] for r in order}}
    for g,rows in per.items():
        record=SeqIO.read(data/f'regions/{g}.gbk','genbank')
        fs={next(f.qualifiers[k][0] for k in ['locus_tag','gene','ID'] if k in f.qualifiers):f for f in record.features if f.type=='CDS' and 'translation' in f.qualifiers}
        for r in rows:
            f=fs[r['target_locus']]
            assert (int(r['start'])-1,int(r['end']),int(r['strand']))==(int(f.location.start),int(f.location.end),f.location.strand)
    return expected,per

def check_artists(fig,artists,bars,axb,scale):
    ref,per=independent_sources()
    selected=sorted([g for g,r in per.items() if len(r)>=43],key=lambda g:(g!='BAA554_hybrid',g))
    source_cons={r['gene']:int(r['genomes_with_ortholog_of_40']) for r in csv.DictReader(io.StringIO((ROOT/'data/region002_gene_conservation.csv').read_text()))}
    assert len(artists)==45+sum(len(per[g]) for g in selected)
    for panel,g,patch in artists:
        points=patch.get_xy();xs=points[:,0];ys=points[:,1]
        if panel=='a':
            lo,hi,st,_=ref[g['gene']];y=0
        else:
            rows=per[g['genome']]
            row=next(r for r in rows if r['target_locus']==g['target_locus'])
            seq=sorted(rows,key=lambda r:int(r['start']))
            reverse=int(seq[-1]['ref_gene_order'])<int(seq[0]['ref_gene_order'])
            start=int(row['start'])-1;end=int(row['end'])
            lo0=min(int(r['start'])-1 for r in rows);hi0=max(int(r['end']) for r in rows)
            lo,hi=(hi0-end,hi0-start) if reverse else (start-lo0,end-lo0)
            st=int(row['strand'])*(-1 if reverse else 1);y=-selected.index(g['genome'])
        assert np.isclose(xs.min(),lo/1000) and np.isclose(xs.max(),hi/1000)
        assert np.isclose(points[2,0],hi/1000 if st==1 else lo/1000) and np.isclose(points[2,1],y)
        assert abs(points[2,0]-points[1,0])>0
    assert len(bars)==45
    kcb=(ROOT/'data/knownclusterblast/BAA554_hybrid.txt').read_text().split('>>')[1]
    tags=set(re.findall(r'^(ctg1_\d+)\t',kcb,re.M))
    expected_ticks=sorted((coords[0]+coords[1])/2000 for gene,coords in ref.items() if gene.replace('BAA554_hybrid_contig_1_','ctg1_') in tags)
    actual_ticks=sorted(float(line.get_xdata()[0]) for line in fig.axes[0].lines)
    assert len(actual_ticks)==20 and np.allclose(actual_ticks,expected_ticks)
    for bar,(gene,coords) in zip(bars,sorted(ref.items(),key=lambda kv:kv[1][3])):
        assert np.isclose(bar.get_x()+bar.get_width()/2,coords[3])
        assert bar.get_height()==source_cons[gene]
    # Display transformations must keep the scale bar in the gene-track units.
    fig.canvas.draw()
    assert np.allclose(axb.transData.transform([(0,0),(10,0)])[:,0],scale.transData.transform([(0,0),(10,0)])[:,0])
    renderer=fig.canvas.get_renderer()
    texts=[t for ax in fig.axes for t in ax.texts]+list(fig.texts)
    texts += [t for ax in fig.axes if ax.axison for t in [*ax.get_xticklabels(),*ax.get_yticklabels(),ax.xaxis.label,ax.yaxis.label]]
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        assert t.get_fontsize()>=8,(t.get_text(),t.get_fontsize())
        box=t.get_window_extent(renderer)
        assert box.x0>=-1 and box.y0>=-1 and box.x1<=fig.bbox.width+1 and box.y1<=fig.bbox.height+1,('Clipped text',t.get_text())
    return selected,per

def raster(path):
    doc=pypdfium2.PdfDocument(str(path));result=np.asarray(doc[0].render(scale=2).to_pil());doc.close();return result

def validate(output,pdffonts):
    fig,audit,artists,bars,axb,scale=build()
    selected,per=check_artists(fig,artists,bars,axb,scale)
    assert json.loads((output/'plot_data.json').read_text())==audit
    with pdfplumber.open(output/'figure4.pdf') as doc:
        assert len(doc.pages)==1;page=doc.pages[0];w,h=page.width,page.height
        assert np.allclose([w,h],np.array([170,220])*72/25.4,atol=.001)
        assert not [c for c in page.chars if c['x0']<-.01 or c['x1']>w+.01 or c['top']<-.01 or c['bottom']>h+.01]
        assert all(c['size']>=7.99 for c in page.chars if c['upright'] and c['text'].strip())
        text=page.extract_text()
        for t in audit['tracks']:assert t['label'] in text,t['label']
        for name in MARK.values():assert name in text
        for i,t in enumerate(audit['tracks']):
            center=h-(.32+.405*((-i)-(-len(selected)+.3))/(.7-(-len(selected)+.3)))*h
            band=page.crop((.812*w,center-5,.985*w,center+5)).extract_text()
            assert f"{len(t['genes'])}/45" in band and f"{t['span_bp']/1000:.1f}" in band,(t['genome'],band)
        scale_y=(1-.305)*h
        lines=[l for l in page.lines if abs(l['top']-scale_y)<.1 and abs(l['bottom']-scale_y)<.1]
        expected=10/86*.505*w
        assert any(abs(l['width']-expected)<.001 for l in lines),'Exported 10 kb scale incorrect'
        # Each exported gene polygon must match the source-derived in-memory polygon.
        curves=[np.array(c['pts']) for c in page.curves if len(c.get('pts',[]))==6]
        for _,_,patch in artists:
            points=patch.get_transform().transform(patch.get_xy())*72/fig.dpi
            points[:,1]=h-points[:,1]
            assert any(c.shape==points.shape and np.allclose(c,points,atol=.001) for c in curves),'Missing/altered exported gene polygon'
    font=subprocess.run([pdffonts,str(output/'figure4.pdf')],capture_output=True,text=True,check=True)
    assert not font.stderr.strip(),font.stderr
    assert 'TrueType' in font.stdout and 'Type 0C' not in font.stdout
    with tempfile.TemporaryDirectory(prefix='figure4-validation-') as tmp:
        p=Path(tmp);fig.savefig(p/'figure4.pdf');fig.savefig(p/'figure4.png',dpi=600)
        assert np.array_equal(raster(p/'figure4.pdf'),raster(output/'figure4.pdf'))
        assert np.array_equal(np.asarray(Image.open(p/'figure4.png')),np.asarray(Image.open(output/'figure4.png')))
    plt.close(fig)
    print(f'PASS: {len(artists)} source-derived gene arrows, 45 conservation bars, {len(selected)} row counts/spans, PDF polygons, 10 kb scale, embedded TrueType fonts, 170 x 220 mm, >=8 pt upright text, and identical fresh PDF/PNG rendering.')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--pdffonts',default='pdffonts');a=p.parse_args();validate(a.output_dir,a.pdffonts)
