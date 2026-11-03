#!/usr/bin/env python3
"""Derive Figure 4c summaries from the selected complete run and verify plot inputs."""
import argparse,collections,csv,hashlib,json,statistics,subprocess,sys
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--module',type=Path,required=True,help='analyses/figure4 directory')
    p.add_argument('--output',type=Path,required=True,help='New or empty audit output directory')
    a=p.parse_args();m=a.module.resolve();out=a.output.resolve()
    if out.exists() and any(out.iterdir()):p.error('Output directory must be new or empty')
    results=m/'reconstruction/full_panel_repeat/results'
    validator=Path(__file__).resolve().parents[1]/'reconstruction/bin/validate_evidence.py'
    raw=subprocess.run([sys.executable,str(validator),'--results',str(results),'--quiet'],
                       capture_output=True,text=True)
    if raw.returncode:
        raise ValueError('Raw evidence validation failed: '+raw.stdout+raw.stderr)
    calls_path=results/'rbh_calls.tsv';json_path=results/'rbh_calls.json'
    json_rows=json.loads(json_path.read_text())
    queries=set(json.loads((results/'query_id_check.json').read_text()))
    genome_ids=set(json.loads((results/'proteome_manifest.json').read_text()))
    if len(queries)!=45 or len(genome_ids)!=40:
        raise ValueError('Expected manifests for 45 queries and 40 genomes')
    expected={(g,q) for g in genome_ids for q in queries}
    with calls_path.open() as f:
        reader=csv.DictReader(f,delimiter='\t');rows=list(reader)
        fields=reader.fieldnames
    keys=[(v['genome'],v['query']) for v in rows]
    if len(keys)!=len(expected) or set(keys)!=expected:
        raise ValueError('TSV must cover exactly the manifest 40 by 45 decision grid')
    canonical={(v['genome'],v['query']):v for v in json_rows}
    if len(json_rows)!=len(expected) or set(canonical)!=expected:
        raise ValueError('JSON must cover exactly the manifest decision grid')
    for row in rows:
        source=canonical[(row['genome'],row['query'])]
        if len(fields)!=len(source) or set(fields)!=set(source):
            raise ValueError('TSV/JSON field names differ')
        # DictWriter serializes None as empty and all other saved values via str.
        if any(row[k]!=('' if v is None else str(v)) for k,v in source.items()):
            raise ValueError('TSV/JSON field mismatch: '+row['genome']+' / '+row['query'])
    if any(r['reciprocal'] not in ('True','False') for r in rows):
        raise ValueError('Invalid decision')
    accepted=[r for r in rows if r['reciprocal']=='True']
    # run_rbh.predict() creates IDs as {contig}_{1-based ordinal}.
    # Bind contig counts to validated target IDs, not an independently editable field.
    for row in accepted:
        contig,separator,ordinal=row['forward_hit'].rpartition('_')
        if not separator or not ordinal.isdigit() or int(ordinal)<1 or row['contig']!=contig:
            raise ValueError('Accepted contig disagrees with validated protein ID')
    genes=collections.Counter(r['query'] for r in accepted)
    genomes=collections.Counter(r['genome'] for r in accepted)
    gene_path=m/'data/region002_gene_conservation.csv';genome_path=m/'data/region002_rbh_matrix.csv'
    with gene_path.open() as f:gene_input=list(csv.DictReader(f))
    with genome_path.open() as f:genome_input=list(csv.DictReader(f))
    if len(gene_input)!=45 or len(genome_input)!=40 or {v['gene'] for v in gene_input}!=queries or {v['genome'] for v in genome_input}!=genome_ids:raise ValueError('Unexpected input dimensions')
    derived_gene=[];derived_genome=[]
    for r in gene_input:
        n=genes[r['gene']]
        if n!=int(r['genomes_with_ortholog_of_40']):raise ValueError('Gene count mismatch: '+r['gene'])
        derived_gene.append(dict(gene_order=r['gene_order'],gene=r['gene'],genomes_with_ortholog_of_40=n))
    contigs={}
    for r in genome_input:
        g=r['genome'];hits=[v for v in accepted if v['genome']==g]
        contigs[g]=len({v['contig'] for v in hits})
        if genomes[g]!=int(r['rbh_of_45']) or contigs[g]!=int(r['n_contigs']):raise ValueError('Genome count mismatch: '+g)
        derived_genome.append(dict(genome=g,rbh_of_45=genomes[g],n_contigs=contigs[g]))
    if set(genes)!={r['gene'] for r in gene_input} or set(genomes)!={r['genome'] for r in genome_input}:raise ValueError('Identifier mismatch')
    out.mkdir(parents=True,exist_ok=True)
    for name,data in [('derived_gene_conservation.csv',derived_gene),('derived_genome_recovery.csv',derived_genome)]:
        with (out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    report={'status':'PASS','source':'reconstruction/full_panel_repeat/results/rbh_calls.tsv','calls':len(rows),'accepted':len(accepted),'genes':len(genes),'genomes':len(genomes),'gene_counts_match':True,'genome_and_contig_counts_match':True,'universal_genes':sum(n==40 for n in genes.values()),'variable_counts':sorted(n for n in genes.values() if n!=40),'per_genome_count_distribution':dict(collections.Counter(genomes.values())),'median_per_genome':statistics.median(genomes.values()),'contig_distribution':dict(collections.Counter(contigs.values())),'input_sha256':{str(path.relative_to(m)):digest(path) for path in (calls_path,json_path,gene_path,genome_path,results/'query_id_check.json',results/'proteome_manifest.json',results/'forward_hits.json',results/'reverse_hits.json',results/'environment.json',results/'presence_matrix.tsv')},'output_sha256':{path.name:digest(path) for path in sorted(out.glob('*.csv'))},'raw_evidence_validated':True,'tsv_json_all_fields_match':True,'manifest_grid_complete':True,'accepted_contigs_match_protein_ids':True,'scope':'Raw-hit validation is required by this command, followed by exact keyed TSV/JSON agreement, complete manifest-grid validation, and comparison with plotted counts.'}
    (out/'pair-to-plot-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
