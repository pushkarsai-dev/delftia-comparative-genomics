"""Validated, explicit provenance for each of the 440 displayed cells."""
import csv
import hashlib
import json
from pathlib import Path
import re

MODULE = Path(__file__).resolve().parents[1]
FAMILIES = [
    ('PF00282', 'Group II PLP-DC', 'Target families'),
    ('PF18631', 'Cucumopine_C', 'Target families'),
    ('PF01977', 'UbiD', 'UbiD/UbiX-related'),
    ('PF02441', 'UbiX/CoaC-like', 'UbiD/UbiX-related'),
    ('PF01276', 'Basic-aa DC', 'Other decarboxylases'),
    ('PF02784', 'DAP/Orn DC', 'Other decarboxylases'),
    ('PF00290', 'Trp synthase alpha', 'Trp metabolism'),
    ('PF00218', 'IGPS', 'Trp metabolism'),
    ('PF03301', 'Trp 2,3-dioxygenase', 'Trp metabolism'),
    ('PF03641', 'LOG family', 'Search controls'),
    ('PF02668', 'TauD/TfdA', 'Search controls'),
]
STUDY = {'ALG19_illumina':'ALG19', 'ATCC13751_hybrid':'ATCC 13751',
         'ATCC15668_hybrid':'ATCC 15668', 'BAA554_hybrid':'BAA-554'}
EXCEPTION = 'GCA_055379645.1'
EXPECTED = {
    'PF00282':(0,0,0), 'PF18631':(0,0,0), 'PF01977':(1,0,1),
    'PF02441':(41,1,2), 'PF01276':(45,1,2), 'PF02784':(40,1,1),
    'PF00290':(40,1,1), 'PF00218':(40,1,1), 'PF03301':(80,2,2),
    'PF03641':(80,2,2), 'PF02668':(202,4,8),
}

def load_data(data_dir=None, verify_hashes=True):
    data_dir = Path(data_dir) if data_dir else MODULE/'data'
    if verify_hashes:
        hashes = json.loads((data_dir/'input_checksums.json').read_text())
        for name, expected in hashes.items():
            if Path(name).name != name:
                raise ValueError('invalid checksum path')
            if hashlib.sha256((data_dir/name).read_bytes()).hexdigest() != expected:
                raise ValueError(f'input checksum mismatch: {name}')
    with (data_dir/'panel40_family_counts.csv').open(newline='') as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        rows = list(reader)
    if len(rows) != 40 or len({r['genome'] for r in rows}) != 40:
        raise ValueError('expected 40 unique genomes')
    if not set(STUDY).issubset({r['genome'] for r in rows}) or EXCEPTION not in {r['genome'] for r in rows}:
        raise ValueError('required genome missing')
    if len(columns) != len(set(columns)):
        raise ValueError('duplicate CSV column')
    for r in rows:
        for key, value in r.items():
            if key != 'genome' and (value is None or not re.fullmatch(r'0|[1-9][0-9]*', value)):
                raise ValueError(f'nonnegative integer required: {key}')
    if sum(int(r['n_proteins']) for r in rows) != 239999:
        raise ValueError('protein total differs')
    # This is an explicitly recorded separate negative result, not a fillna(0).
    evidence = (data_dir/'pictet_spenglerase_search.txt').read_text()
    if not re.search(r'PF18631\s+Cucumopine_C[^\n]*0 hits in 40 of 40 genomes', evidence):
        raise ValueError('separate PF18631 negative-result evidence missing')
    rows.sort(key=lambda r:(r['genome'] not in STUDY,r['genome']))
    keys = {}
    for fam,_,_ in FAMILIES:
        if fam == 'PF18631':
            continue
        matches = [c for c in columns if c.startswith(fam+'_')]
        if len(matches) != 1:
            raise ValueError(f'exactly one column required: {fam}')
        keys[fam] = matches[0]
    values = [[0 if fam == 'PF18631' else int(r[keys[fam]]) for fam,_,_ in FAMILIES] for r in rows]
    independent = json.loads((data_dir/'independent_search.json').read_text())
    if independent['n_proteins'] != 239999 or set(independent['counts']['per_genome']) != {r['genome'] for r in rows}:
        raise ValueError('independent search panel differs')
    for r, vals in zip(rows,values):
        fresh = independent['counts']['per_genome'][r['genome']]
        for j,(fam,_,_) in enumerate(FAMILIES):
            if type(fresh[fam]) is not int or fresh[fam] != vals[j]:
                raise ValueError(f'independent count disagreement: {r["genome"]} {fam}')
        # Use the explicitly recorded fresh result for this separately searched family.
        vals[1] = fresh['PF18631']
    for j,(fam,_,_) in enumerate(FAMILIES):
        vals = [r[j] for r in values]
        if (sum(vals),min(vals),max(vals)) != EXPECTED[fam]:
            raise ValueError(f'archived family summary differs: {fam}')
    for r,vals in zip(rows,values):
        if bool(vals[2]) != (r['genome'] == EXCEPTION) or (vals[3] == 2) != (r['genome'] == EXCEPTION):
            raise ValueError('UbiD/PF02441 exception moved')
    return [r['genome'] for r in rows], values

def reserve_output(path):
    path = Path(path).resolve()
    for protected in [MODULE/'data',MODULE/'scripts',MODULE/'tests',MODULE/'docs',MODULE/'reference']:
        protected = protected.resolve()
        if path == protected or path in protected.parents or protected in path.parents:
            raise ValueError('output overlaps protected module content')
    if path.exists():
        raise ValueError('output directory already exists; choose a new path')
    path.mkdir(parents=True,exist_ok=False)
    return path
