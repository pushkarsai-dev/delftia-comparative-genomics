"""Repeat region-restricted reciprocal searches; this is not the whole-proteome run."""
import argparse, csv, hashlib, json, re
from pathlib import Path
import pyhmmer
from Bio import SeqIO

ROOT=Path(__file__).resolve().parents[1]
EV=1e-5

def read_cds(path, mode='feature'):
    if mode=='historical':
        out=[]
        # Archived extraction split at CDS lines but searched the entire block,
        # including following domain features, for the first locus tag.
        for b in path.read_text().split('\n     CDS             ')[1:]:
            loc=b.split('\n',1)[0].strip()
            coords=re.findall(r'(\d+)\.\.>?(\d+)',loc)
            tag=re.search(r'/locus_tag="([^"]+)"',b)
            tr=re.search(r'/translation="([^"]+)"',b,re.S)
            if coords and tag and tr:
                out.append((tag[1],int(coords[0][0]),int(coords[-1][1]),-1 if loc.startswith('complement') else 1,re.sub(r'\s+','',tr[1])))
        return out
    rec=SeqIO.read(path,'genbank');out=[]
    for f in rec.features:
        if f.type!='CDS' or 'translation' not in f.qualifiers:continue
        tag=next((f.qualifiers[k][0] for k in ['locus_tag','gene','ID'] if k in f.qualifiers),None)
        if tag is None:raise ValueError('CDS lacks an identifier')
        out.append((tag,int(f.location.start)+1,int(f.location.end),f.location.strand,f.qualifiers['translation'][0]))
    return out

def best(queries,targets):
    alpha=pyhmmer.easel.Alphabet.amino()
    def block(rows):
        return pyhmmer.easel.DigitalSequenceBlock(alpha,[pyhmmer.easel.TextSequence(name=r[0].encode(),sequence=r[4]).digitize(alpha) for r in rows])
    out={}
    for hits in pyhmmer.phmmer(block(queries),block(targets),E=EV,cpus=1):
        candidates=sorted((h for h in hits if h.evalue<EV),key=lambda h:h.evalue)
        if candidates:
            q=hits.query.name;h=candidates[0]
            dec=lambda s:s.decode() if isinstance(s,bytes) else s
            out[dec(q)]=(dec(h.name),h.evalue,h.score)
    return out

def run(output,mode):
    if output.exists():raise ValueError('Output directory must not already exist')
    from figure4_data import verify_inputs
    verify_inputs(ROOT/'data')
    manifest=json.loads((ROOT/'data/region_manifest.json').read_text())
    ref=sorted(read_cds(ROOT/'data/regions/BAA554_hybrid.gbk',mode),key=lambda r:r[1])
    assert len(ref)==45
    output.mkdir(parents=True)
    rows=[];audit=[]
    for item in manifest:
        genome=item['genome'];path=ROOT/'data'/item['genbank']
        target=read_cds(path,mode)
        fwd=best(ref,target);back=best(target,ref)
        info={r[0]:r for r in target};n=0
        for order,q in enumerate(ref,1):
            hit=fwd.get(q[0])
            if hit is None or back.get(hit[0],(None,))[0]!=q[0]:continue
            t,s,e,st,_=info[hit[0]];n+=1
            rows.append(dict(genome=genome,record=item['record'],ref_gene_order=order,ref_locus=q[0],target_locus=t,start=s,end=e,strand=st))
        audit.append(dict(genome=genome,genbank_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),targets=len(target),rbh=n,forward=fwd,reverse=back))
        (output/'search_progress.json').write_text(json.dumps(audit,indent=2)+'\n')
        print(f'{genome}: {n}/45 reciprocal matches ({len(target)} target CDS)',flush=True)
    with (output/'region_gene_order.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    with (ROOT/'data/region_gene_order.tsv').open() as handle:
        old=list(csv.DictReader(handle,delimiter='\t'))
    def rowkey(r):return tuple(str(r[k]) for k in ['genome','record','ref_gene_order','ref_locus','target_locus','start','end','strand'])
    before={rowkey(r) for r in old};after={rowkey(r) for r in rows}
    summary=dict(pyhmmer=pyhmmer.__version__,E=EV,parser_mode=mode,genomes=len(audit),rows=len(rows),
                 archived_rows=len(old),missing=sorted(before-after),added=sorted(after-before),
                 limitation='Does not repeat panel-wide whole-proteome RBH search used for conservation counts.')
    (output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['missing','added']},indent=2),flush=True)
    print(f'Missing archived rows: {len(before-after)}; additional rows: {len(after-before)}',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',required=True,type=Path);p.add_argument('--parser',choices=['historical','feature'],default='feature');a=p.parse_args()
    run(a.output_dir,a.parser)
