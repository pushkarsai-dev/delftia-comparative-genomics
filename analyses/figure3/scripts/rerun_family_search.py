"""Independent, read-only-input reproduction of the archived family counts.

New outputs only. No historical cache is loaded. Both pooled and per-genome
target databases are evaluated because the original summary does not identify
which target arrangement generated the family-count CSV.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import pyhmmer
import pyrodigal

FAMILIES = ['PF00282','PF01977','PF02441','PF03088','PF20067','PF01276',
            'PF02784','PF02668','PF03641','PF03301','PF00155','PF00290',
            'PF00218','PF18631','PF02317','PF10100','PF02423']

def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for b in iter(lambda: handle.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()

def fasta(path):
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)

def string(value):
    return value.decode() if isinstance(value, bytes) else str(value)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project', type=Path, required=True)
    ap.add_argument('--pfam', type=Path, required=True)
    ap.add_argument('--outdir', type=Path, required=True)
    args = ap.parse_args()
    if args.outdir.exists():
        ap.error('output directory must not exist')
    if pyrodigal.__version__ != '3.7.1' or pyhmmer.__version__ != '0.12.2':
        ap.error('requires Pyrodigal 3.7.1 and PyHMMER 0.12.2')
    panel = args.project / 'A3 antiSMASH Local (2 Aug 2026)/panel40_fasta'
    a10 = args.project / 'A10 InterProScan (14 Aug 2026)'
    files = sorted(panel.glob('*.fna'))
    archived = list(csv.DictReader((a10/'panel40_family_counts.csv').open()))
    archive = {r['genome']:r for r in archived}
    recorded = json.loads((a10/'Deposition package (8 Sep 2026)/reference_run/results/01_proteomes_manifest.json').read_text())
    assert len(files) == len(archive) == 40
    assert {p.stem for p in files} == set(archive)
    input_hashes = {p.name:digest(p) for p in files}
    assert input_hashes == recorded['input_sha256'], 'input FASTAs changed'
    pfam_md5 = digest(args.pfam, 'md5')
    assert pfam_md5 == '2cacd81054f0c267718e752ff9322bf1', 'Pfam database changed'
    args.outdir.mkdir(parents=True)
    profiles = {}
    with pyhmmer.plan7.HMMFile(str(args.pfam)) as handle:
        for h in handle:
            acc = string(h.accession).split('.')[0]
            if acc in FAMILIES:
                profiles[acc] = h
    assert set(profiles) == set(FAMILIES)
    profile_dir = args.outdir/'profiles'
    profile_dir.mkdir()
    for acc, h in profiles.items():
        with (profile_dir/f'{acc}.hmm').open('wb') as out:
            h.write(out)
    abc = pyhmmer.easel.Alphabet.amino()
    proteins = []
    per_genome = {}
    start = time.monotonic()
    with (args.outdir/'predicted_proteins.faa').open('w') as faa:
        for path in files:
            seqs = [s.encode() for _,s in fasta(path)]
            finder = pyrodigal.GeneFinder(meta=False)
            finder.train(*seqs)
            block = pyhmmer.easel.DigitalSequenceBlock(abc)
            for i, seq in enumerate(seqs):
                for j, gene in enumerate(finder.find_genes(seq)):
                    name = f'{path.stem}|{i}_{j}'
                    protein = gene.translate()
                    block.append(pyhmmer.easel.TextSequence(name=name.encode(), sequence=protein).digitize(abc))
                    faa.write(f'>{name}\n{protein}\n')
            per_genome[path.stem] = block
            proteins.extend(block)
            assert len(block) == int(archive[path.stem]['n_proteins'])
            print(path.stem, len(block), flush=True)
    assert len(proteins) == 239999
    print('Prediction seconds', time.monotonic()-start, flush=True)
    all_counts, hit_records = {}, []
    for mode in ['per_genome', 'pooled']:
        counts = {g:{f:0 for f in FAMILIES} for g in archive}
        databases = per_genome.items() if mode == 'per_genome' else [('panel', pyhmmer.easel.DigitalSequenceBlock(abc, proteins))]
        for dbname, block in databases:
            for hits in pyhmmer.hmmer.hmmsearch(list(profiles.values()), block, E=1e-3, cpus=4):
                family = string(hits.query.accession).split('.')[0]
                for hit in hits:
                    if hit.evalue < 1e-3:
                        target = string(hit.name)
                        genome = target.split('|')[0]
                        counts[genome][family] += 1
                        hit_records.append(dict(mode=mode, target_db=dbname, family=family,
                                                target=target, evalue=float(hit.evalue), score=float(hit.score)))
        all_counts[mode] = counts
        disagreements = []
        for g,row in archive.items():
            for key,value in row.items():
                if key.startswith('PF') and counts[g][key.split('_')[0]] != int(value):
                    disagreements.append(dict(genome=g, family=key.split('_')[0], archived=int(value), reproduced=counts[g][key.split('_')[0]]))
        (args.outdir/f'{mode}_comparison.json').write_text(json.dumps(disagreements, indent=2)+'\n')
        print(mode, 'archived-cell disagreements', len(disagreements), flush=True)
    with (args.outdir/'hits.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(hit_records[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(hit_records)
    result = dict(versions=dict(python=sys.version, platform=platform.platform(), pyrodigal=pyrodigal.__version__, pyhmmer=pyhmmer.__version__),
                  input_sha256=input_hashes, pfam_md5=pfam_md5,
                  profiles={f:digest(profile_dir/f'{f}.hmm') for f in FAMILIES},
                  n_proteins=len(proteins), n_per_genome={g:len(b) for g,b in per_genome.items()},
                  counts=all_counts, scope='New reproduction, not recovered original execution code',
                  seconds=time.monotonic()-start)
    (args.outdir/'rerun.json').write_text(json.dumps(result, indent=2)+'\n')
    assert input_hashes == {p.name:digest(p) for p in files}
    print('COMPLETE', flush=True)

if __name__ == '__main__':
    main()
