"""Independently compare exported geometry/text with the archived JSON and TSV."""
import argparse, csv, hashlib, json, math
from pathlib import Path
import pdfplumber
from pdfminer.pdftypes import resolve1
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
ORDER = ('GENE3D', 'SUPERFAMILY', 'PFAM', 'PANTHER')
COLORS = ((196/255,218/255,235/255),(214/255,229/255,239/255),
          (152/255,191/255,219/255),(176/255,203/255,224/255))

def require(condition, message):
    if not condition:
        raise ValueError(message)

def close(a, b, tolerance=0.002):
    return abs(a-b) <= tolerance

def color_is(a, b):
    return isinstance(a,(tuple,list)) and len(a)==len(b) and all(close(x,y,1e-8) for x,y in zip(a,b))

def sources(data):
    """No imports from the plotting implementation or its source parser."""
    manifest=json.loads((data/'input_checksums.json').read_text())
    require(set(manifest)=={'InterProScan_results.json','InterProScan_results.tsv','BAA554_two_loci.faa'},'Input manifest filenames')
    for name,digest in manifest.items():
        require(hashlib.sha256((data/name).read_bytes()).hexdigest()==digest,'Input checksum: '+name)
    seqs={}
    for line in (data/'BAA554_two_loci.faa').read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            require(key not in seqs,'Duplicate FASTA ID');seqs[key]=''
        elif line: seqs[key]+=line
    expected_ids=['HFJBPOEG_02821','HFJBPOEG_02843']
    require(set(seqs)==set(expected_ids),'Protein IDs')
    raw=json.loads((data/'InterProScan_results.json').read_text())
    require(raw['interproscan-version']=='5.78-109.0' and len(raw['results'])==2,'Archive identity')
    with (data/'InterProScan_results.tsv').open() as stream:
        tsv=list(csv.reader(stream,delimiter='\t'))
    require(len(tsv)==8,'TSV record count')
    result=[]
    for pid,length in zip(expected_ids,[289,329]):
        seq=seqs[pid]; require(len(seq)==length,'Protein length')
        candidates=[p for p in raw['results'] if p['sequence']==seq]
        require(len(candidates)==1,'JSON sequence identity');p=candidates[0]
        require(hashlib.md5(seq.encode()).hexdigest()==p['md5'],'Sequence MD5')
        require(len(p['matches'])==4,'Signature count')
        hits=[]
        for library in ORDER:
            matches=[m for m in p['matches'] if m['signature']['signatureLibraryRelease']['library']==library]
            require(len(matches)==1,'Library uniqueness');m=matches[0];sig=m['signature']
            require(len(m['locations'])==1,'Location count');loc=m['locations'][0]
            fragments=sorted((f['start'],f['end']) for f in loc['location-fragments'])
            e=loc['evalue'] if library in ('GENE3D','PFAM') else m['evalue']
            require(type(e) in (int,float) and math.isfinite(e) and e>=0,'E-value')
            rows=[r for r in tsv if r[0]==pid and r[4]==sig['accession']]
            require(len(rows)==1,'TSV identity');r=rows[0]
            require((int(r[6]),int(r[7]),float(r[8]))==(loc['start'],loc['end'],e),'TSV/JSON coordinates or E-value')
            hits.append((sig['accession'],fragments,e))
        result.append((pid,length,hits))
    require(seqs[expected_ids[0]].count('PGGFGTLDE')==1 and seqs[expected_ids[0]].index('PGGFGTLDE')==210,'Motif 211-219')
    require(sum(len(h[1]) for p in result for h in p[2])==9,'Nine matched segments')
    return result

def inspect_pdf(path, expected):
    with pdfplumber.open(path) as pdf:
        require(len(pdf.pages)==1,'One page');p=pdf.pages[0]
        w,h=p.width,p.height
        require(close(w,170/25.4*72) and close(h,160/25.4*72),'Page dimensions 170 x 160 mm')
        chars=p.chars
        require(chars and min(c['size'] for c in chars)>=7.49,'Minimum 7.5-point text')
        require(all(c['x0']>=0 and c['x1']<=w and c['top']>=0 and c['bottom']<=h for c in chars),'Text clipping')
        fonts=resolve1(resolve1(p.page_obj.resources)['Font'])
        for item in fonts.values():
            font=resolve1(item)
            for child in resolve1(font.get('DescendantFonts',[item])):
                child=resolve1(child);desc=resolve1(child.get('FontDescriptor',{}))
                require('FontFile2' in desc,'Font must embed actual TrueType')
                require('FontFile3' not in desc,'No CFF font mismatch')
        text=p.extract_text()
        require(text.count('PGGFGTLDE')==1,'One motif label')
        require(text.count('Residue position')==2,'Axis labels')
        colored=[r for r in p.rects if any(color_is(r['non_stroking_color'],c) for c in COLORS)]
        require(len(colored)==9,'Nine rendered signature segments')
        expected_rects=[]
        x=lambda pos: w*(.205+.565*(pos-.5)/329)
        for i,(pid,length,hits) in enumerate(expected):
            bottom=.565 if i==0 else .085
            y=lambda pos: h*(bottom+.33*(pos+.65)/5.75)
            require(f'{pid} · {length} aa' in text,'Protein title')
            gray=[r for r in p.rects if isinstance(r['non_stroking_color'],(int,float)) and close(r['non_stroking_color'],217/255,1e-8) and close(r['y0'],y(4.26))]
            require(len(gray)==1 and close(gray[0]['x0'],x(.5)) and close(gray[0]['x1'],x(length+.5)),'Full protein geometry')
            for k,(accession,fragments,e) in enumerate(hits):
                require(text.count(accession)==1,'Signature label '+accession)
                rowtext=p.crop((0,h-y(3-k+.80),w,h-y(3-k-.18))).extract_text()
                require(accession in rowtext,'Signature is on the wrong row: '+accession)
                spans=', '.join(f'{a}–{b}' for a,b in fragments)
                require(spans in rowtext,'Displayed coordinates '+accession)
                mant,exponent=f'{e:.1e}'.split('e')
                label='E = '+mant+' × 10'+str(int(exponent)).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
                require(label in rowtext,'Displayed E-value '+accession)
                for a,b in fragments: expected_rects.append((x(a-.5),y(3-k),x(b+.5),y(3-k+.56),COLORS[k]))
        for x0,y0,x1,y1,color in expected_rects:
            matches=[r for r in colored if color_is(r['non_stroking_color'],color) and all(close(r[k],v) for k,v in zip(('x0','y0','x1','y1'),(x0,y0,x1,y1)))]
            require(len(matches)==1,'Segment geometry differs from archived fragments')
        red=[l for l in p.lines if color_is(l['stroking_color'],(178/255,34/255,34/255))]
        require(len(red)==1 and close(red[0]['x0'],x(210.5)) and close(red[0]['x1'],x(219.5)),'Motif geometry 211-219')
        connectors=[l for l in p.lines if color_is(l['stroking_color'],(82/255,105/255,123/255))]
        require(len(connectors)==1 and close(connectors[0]['x0'],x(40.5)) and close(connectors[0]['x1'],x(282.5)),'Discontinuous match connector')
        return {'signatures':8,'matched_segments':9,'protein_lengths':[289,329], 'motif':[211,219],'minimum_font_pt':min(c['size'] for c in chars),'page_mm':[170,160],'embedded_truetype':True}

def validate(output, data=ROOT/'data'):
    expected=sources(data)
    report=inspect_pdf(output/'figure5.pdf',expected)
    plotted=json.loads((output/'plotted_data.json').read_text())
    require(len(plotted)==2,'Plotted record count')
    for p,(pid,length,hits) in zip(plotted,expected):
        require(p['id']==pid and p['length']==length and len(p['hits'])==4,'Plotted protein identity')
        for row,(accession,fragments,e) in zip(p['hits'],hits):
            require(row['accession']==accession and row['fragments']==[list(f) for f in fragments] and row['evalue']==e,'Plotted data mismatch')
    require(plotted[0]['motif']=={'sequence':'PGGFGTLDE','start':211,'end':219},'Plotted motif')
    with Image.open(output/'figure5.png') as png:
        require(png.size==(4015,3779),'600 dpi PNG dimensions')
    caption=(ROOT/'docs/caption-and-alt-text.md').read_text()
    alt=caption.split('## Alt text\n\n',1)[1].split('\n\n',1)[0]
    require(len(alt)<=255 and '\n' not in alt,'Alt text maximum 255 characters')
    for phrase in ('41-84 and 115-282','not biochemical validation','211-219','common residue scale'):
        require(phrase in caption,'Caption lacks required detail: '+phrase)
    report['alt_text_characters']=len(alt)
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output_dir',type=Path);ap.add_argument('--data-dir',type=Path,default=ROOT/'data')
    args=ap.parse_args()
    try: print(json.dumps(validate(args.output_dir,args.data_dir),indent=2))
    except (ValueError,KeyError,TypeError,OSError) as e: ap.exit(2,'FAILED: '+str(e)+'\n')
