"""Validated inputs and coordinate transforms for Figure 4."""
import csv, hashlib, json, re, io
from pathlib import Path
from collections import defaultdict

ROOT=Path(__file__).resolve().parents[1]
MARK={'BAA554_hybrid_contig_1_1928':'HFJBPOEG_02821','BAA554_hybrid_contig_1_1906':'HFJBPOEG_02843'}
PRETTY={'BAA554_hybrid':'BAA-554','ATCC13751_hybrid':'ATCC 13751'}

def verify_inputs(data):
    manifest=json.loads((data/'input_checksums.json').read_text())
    for name,digest in manifest.items():
        path=data/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Input checksum mismatch: {name}')

def normalize_track(rows):
    rows=sorted(rows,key=lambda r:r['start'])
    orders=[r['ref_gene_order'] for r in rows]
    reverse=orders[-1]<orders[0]
    if (orders[::-1] if reverse else orders)!=sorted(orders):raise ValueError('Non-collinear gene order')
    lo=min(r['start'] for r in rows)-1;hi=max(r['end'] for r in rows)
    tracks=[]
    for r in rows:
        s=r['start']-1;e=r['end']
        x0,x1=(hi-e,hi-s) if reverse else (s-lo,e-lo)
        tracks.append(dict(r,x0=x0,x1=x1,display_strand=r['strand']*(-1 if reverse else 1)))
    return tracks,reverse,hi-lo

def arrow_vertices(x0,x1,y,strand,height=.52,head_limit=.9):
    if not x1>x0 or strand not in [-1,1]:raise ValueError('Invalid arrow interval or strand')
    head=min((x1-x0)*.45,head_limit)
    if strand==1:return [(x0,y-height/2),(x1-head,y-height/2),(x1,y),(x1-head,y+height/2),(x0,y+height/2)]
    return [(x1,y-height/2),(x0+head,y-height/2),(x0,y),(x0+head,y+height/2),(x1,y+height/2)]

def load_data(data=ROOT/'data',verify=True):
    if verify:verify_inputs(data)
    conservation=list(csv.DictReader(io.StringIO((data/'region002_gene_conservation.csv').read_text())))
    cons={r['gene']:int(r['genomes_with_ortholog_of_40']) for r in conservation}
    if len(cons)!=45 or len(conservation)!=45:raise ValueError('Expected 45 unique reference genes')
    genes=[]
    for line in (data/'region002.faa').read_text().splitlines():
        if not line.startswith('>'):continue
        m=re.fullmatch(r'>(\S+)\s+(\d+)-(\d+)\s+strand=(-?1)',line)
        if not m:raise ValueError('Malformed FASTA coordinates')
        genes.append(dict(gene=m[1],start=int(m[2]),end=int(m[3]),strand=int(m[4])))
    genes.sort(key=lambda g:g['start'])
    if len(genes)!=45 or {g['gene'] for g in genes}!=set(cons):raise ValueError('Reference identifier mismatch')
    kcbtext=(data/'knownclusterblast/BAA554_hybrid.txt').read_text()
    blocks=kcbtext.split('>>')[1:]
    if not blocks or not re.match(r'\s*1\. BGC0000984\.',blocks[0]):raise ValueError('Unexpected top cluster')
    matches={line.split('\t')[0] for line in blocks[0].splitlines() if line.startswith('ctg1_')}
    if len(matches)!=20:raise ValueError('Expected 20 unique reference-cluster matches')
    for i,g in enumerate(genes,1):
        g.update(order=i,conservation=cons[g['gene']],kcb=g['gene'].replace('BAA554_hybrid_contig_1_','ctg1_') in matches)
        if not 0<=g['conservation']<=40:raise ValueError('Invalid conservation count')
    per=defaultdict(list)
    rows=list(csv.DictReader(io.StringIO((data/'region_gene_order_corrected.tsv').read_text()),delimiter='\t'))
    seen=set()
    for raw in rows:
        r=dict(raw)
        for k in ['start','end','strand','ref_gene_order']:r[k]=int(r[k])
        key=(r['genome'],r['ref_gene_order'])
        if key in seen:raise ValueError('Duplicate genome/reference pair')
        seen.add(key)
        if not 1<=r['ref_gene_order']<=45 or r['start']<1 or r['end']<r['start'] or r['strand'] not in [-1,1]:raise ValueError('Invalid gene coordinate')
        expected=genes[r['ref_gene_order']-1]['gene'].replace('BAA554_hybrid_contig_1_','ctg1_')
        if r['ref_locus']!=expected:raise ValueError('Reference-order mismatch')
        per[r['genome']].append(r)
    if len(per)!=40:raise ValueError('Expected 40 genomes')
    selection=sorted([g for g,v in per.items() if len(v)>=43],key=lambda g:(g!='BAA554_hybrid',g))
    if len(selection)!=20:raise ValueError('Expected 20 selected genomes')
    tracks=[]
    for g in selection:
        if len({r['record'] for r in per[g]})!=1:raise ValueError('Selected genome spans multiple records')
        normalized,reverse,span=normalize_track(per[g])
        for r in normalized:
            if r['display_strand']!=genes[r['ref_gene_order']-1]['strand']:raise ValueError('Strand disagreement')
        tracks.append(dict(genome=g,label=PRETTY.get(g,g),reversed=reverse,span_bp=span,genes=normalized))
    totals=list(csv.DictReader(io.StringIO((data/'region002_rbh_matrix.csv').read_text())))
    if len(totals)!=40 or {r['genome'] for r in totals}!=set(per):raise ValueError('Panel identifiers differ')
    if sum(int(r['rbh_of_45']) for r in totals)!=sum(cons.values()):raise ValueError('Conservation marginals disagree')
    return genes,tracks,totals
