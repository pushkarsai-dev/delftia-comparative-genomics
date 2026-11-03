"""Validate and reconcile the two archived InterProScan representations."""
from pathlib import Path
import csv,hashlib,json,math

ROOT=Path(__file__).resolve().parents[1]
IDS=('HFJBPOEG_02821','HFJBPOEG_02843')
ORDER=('Gene3D','SUPERFAMILY','Pfam','PANTHER')
LIBRARIES={'GENE3D':('Gene3D','4.3.0'),'SUPERFAMILY':('SUPERFAMILY','1.75'),'PFAM':('Pfam','38.2'),'PANTHER':('PANTHER','19.0')}

def number(x):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0:
        raise ValueError('Invalid nonnegative finite E-value')
    return x

def interval(a,b,length):
    if type(a)!=int or type(b)!=int or not 1<=a<=b<=length:
        raise ValueError('Invalid 1-based inclusive residue interval')
    return [a,b]

def load_data(data_dir=None,verify_hashes=True):
    data_dir=Path(data_dir) if data_dir is not None else ROOT/'data'
    if verify_hashes:
        hashes=json.loads((data_dir/'input_checksums.json').read_text())
        if set(hashes)!={'BAA554_two_loci.faa','InterProScan_results.tsv','InterProScan_results.json'}:
            raise ValueError('Unexpected input manifest')
        for name,h in hashes.items():
            if hashlib.sha256((data_dir/name).read_bytes()).hexdigest()!=h:raise ValueError('Input checksum mismatch: '+name)
    seqs={};key=None
    for line in (data_dir/'BAA554_two_loci.faa').read_text().splitlines():
        if line.startswith('>'):
            key=line[1:].split()[0]
            if key in seqs:raise ValueError('Duplicate FASTA identifier')
            seqs[key]=''
        elif line:
            if key is None or any(c not in 'ACDEFGHIKLMNPQRSTVWY' for c in line):raise ValueError('Invalid protein FASTA')
            seqs[key]+=line
    if set(seqs)!=set(IDS) or len(set(seqs.values()))!=2:raise ValueError('Unexpected proteins')
    raw=json.loads((data_dir/'InterProScan_results.json').read_text())
    if raw['interproscan-version']!='5.78-109.0' or len(raw['results'])!=2:raise ValueError('Unexpected archived version/results')
    json_rows={};seen=set();proteins={}
    for result in raw['results']:
        ids=[p for p,s in seqs.items() if s==result['sequence']]
        if len(ids)!=1 or ids[0] in seen:raise ValueError('JSON sequence missing, duplicate, or different from FASTA')
        pid=ids[0];seen.add(pid);seq=seqs[pid]
        if hashlib.md5(seq.encode()).hexdigest()!=result['md5']:raise ValueError('JSON MD5 mismatch')
        if len(result['matches'])!=4:raise ValueError('Expected four signatures per protein')
        proteins[pid]={'id':pid,'length':len(seq),'sequence_md5':result['md5'],'hits':[]}
        for m in result['matches']:
            sig=m['signature'];release=sig['signatureLibraryRelease']
            db,version=LIBRARIES[release['library']]
            if release['version']!=version or len(m['locations'])!=1:raise ValueError('Unexpected library or locations')
            loc=m['locations'][0];start,end=interval(loc['start'],loc['end'],len(seq))
            frags=sorted(interval(f['start'],f['end'],len(seq)) for f in loc['location-fragments'])
            if not frags or frags[0][0]!=start or frags[-1][1]!=end or any(a[1]>=b[0] for a,b in zip(frags,frags[1:])):raise ValueError('Invalid or overlapping fragments')
            scope='location' if db in ('Gene3D','Pfam') else 'match'
            e=number(loc['evalue'] if scope=='location' else m['evalue'])
            entry=sig.get('entry')
            item={'db':db,'accession':sig['accession'],'start':start,'end':end,'fragments':frags,'evalue':e,'evalue_scope':scope,'interpro':entry['accession'] if entry else None,'description':entry['description'] if entry else sig.get('name'),'library_version':version}
            key=(pid,db)
            if key in json_rows:raise ValueError('Duplicate JSON library')
            json_rows[key]=item
    seen_tsv=set()
    with (data_dir/'InterProScan_results.tsv').open() as f:
        for row in csv.reader(f,delimiter='\t'):
            if len(row)<13:raise ValueError('Truncated TSV row')
            pid,db=row[0],row[3];key=(pid,db)
            if key in seen_tsv or key not in json_rows:raise ValueError('Duplicate or unexpected TSV row')
            seen_tsv.add(key);j=json_rows[key]
            expected=(j['accession'],j['start'],j['end'],j['evalue'],j['interpro'])
            actual=(row[4],int(row[6]),int(row[7]),number(float(row[8])),None if row[11]=='-' else row[11])
            if actual!=expected or int(row[2])!=len(seqs[pid]) or row[1]!=proteins[pid]['sequence_md5'] or row[9]!='T' or row[10]!='14-08-2026':raise ValueError('TSV and JSON disagree: '+str(key))
    if seen_tsv!=set(json_rows):raise ValueError('Missing TSV rows')
    for pid in IDS:
        proteins[pid]['hits']=[json_rows[(pid,db)] for db in ORDER]
    motif='PGGFGTLDE';seq=seqs[IDS[0]]
    if seq.count(motif)!=1:raise ValueError('LOG motif is missing or duplicated')
    start=seq.index(motif)+1
    proteins[IDS[0]]['motif']={'sequence':motif,'start':start,'end':start+len(motif)-1}
    return [proteins[p] for p in IDS]

def e_label(value):
    mantissa,exponent=f'{value:.1e}'.split('e')
    trans=str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    return f'E = {mantissa} × 10{str(int(exponent)).translate(trans)}'
