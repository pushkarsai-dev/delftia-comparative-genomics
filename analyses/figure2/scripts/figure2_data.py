"""Offline recomputation of Figure 2 counts from the archived pangenome object."""
from pathlib import Path
from collections import Counter
from statistics import mean
import csv
import gzip
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = {
    'pangenome.json.gz', 'genome_name_map.tsv', 'families.tsv',
    'presence_absence.tsv', 'node_table.tsv', 'accumulation_curve.tsv',
    'tree_40tips_historical.json', 'tree_40tips_corrected.newick',
    'pan_phylo_table.tsv',
}

def require(condition, message):
    if not condition:
        raise ValueError(message)

def read_json(text):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    def invalid(value):
        raise ValueError(f'Nonfinite JSON value: {value}')
    return json.loads(text, object_pairs_hook=pairs, parse_constant=invalid)

def rows(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames is not None, f'No header: {path.name}')
        require(len(reader.fieldnames) == len(set(reader.fieldnames)), 'Duplicate TSV header')
        result = list(reader)
        require(all(None not in r and all(v is not None for v in r.values()) for r in result),
                f'Ragged TSV: {path.name}')
        return result

def integer(value):
    require(isinstance(value, str) and re.fullmatch(r'\d+', value) is not None,
            f'Not a nonnegative integer: {value!r}')
    return int(value)

def category(n):
    require(type(n) is int and 1 <= n <= 40, 'Occupancy outside 1..40')
    return 'core' if n == 40 else 'soft_core' if n >= 38 else 'shell' if n >= 6 else 'cloud'

def verify_hashes(data):
    data = Path(data)
    pin = read_json((data / 'input_checksums.json').read_text())
    require(isinstance(pin, dict) and set(pin) == REQUIRED, 'Input manifest membership mismatch')
    for name, digest in pin.items():
        require(isinstance(digest, str) and re.fullmatch('[0-9a-f]{64}', digest), 'Invalid digest')
        require((data / name).is_file(), f'Missing input: {name}')
        require(hashlib.sha256((data / name).read_bytes()).hexdigest() == digest,
                f'Checksum mismatch: {name}')
    return pin

def parse_families(obj, names):
    require(isinstance(obj, dict), 'Pangenome must be an object')
    refs = obj['genome_refs']
    require(isinstance(refs, list) and len(refs) == 40 and len(set(refs)) == 40, 'Expected 40 unique genome refs')
    require(set(refs) == set(names) and len(set(names.values())) == 40, 'Genome-map bijection failed')
    result, seen = {}, set()
    for family in obj['orthologs']:
        key = family['id']
        require(isinstance(key, str) and key and key not in result, 'Duplicate/invalid family ID')
        counts = Counter()
        for member in family['orthologs']:
            require(isinstance(member, list) and len(member) == 3, 'Malformed gene membership')
            gene, order, ref = member
            require(isinstance(gene, str) and gene and ref in names, 'Invalid gene or genome ref')
            require(type(order) is int and order >= 0, 'Invalid gene order')
            require((gene, ref) not in seen, 'Gene assigned more than once')
            seen.add((gene, ref))
            counts[names[ref]] += 1
        require(counts, 'Empty family')
        result[key] = counts
    return result, len(seen)

def tree_clades(text):
    """Read the restricted unquoted Newick syntax used by the two pinned trees.

    Preserve polytomies; branch lengths and support do not enter this calculation.
    """
    tokens = re.findall(r'[^\s(),:;]+|[(),:;]', text)
    i, result, seen = 0, [], set()
    def node():
        nonlocal i
        require(i < len(tokens), 'Truncated Newick')
        if tokens[i] == '(':
            i += 1
            leaves = set(node())
            while i < len(tokens) and tokens[i] == ',':
                i += 1
                leaves.update(node())
            require(i < len(tokens) and tokens[i] == ')', 'Unclosed Newick clade')
            i += 1
            result.append(frozenset(leaves))
            if i < len(tokens) and tokens[i] not in ',):;':
                i += 1
        else:
            require(tokens[i] not in ',):;', 'Invalid Newick leaf')
            require(tokens[i] not in seen, 'Duplicate Newick leaf')
            leaves = {tokens[i]}
            seen.update(leaves)
            i += 1
        if i < len(tokens) and tokens[i] == ':':
            require(i + 1 < len(tokens), 'Missing branch length')
            float(tokens[i + 1])
            i += 2
        return leaves
    tips = node()
    require(i == len(tokens) - 1 and tokens[i] == ';', 'Invalid Newick ending')
    return tips, result

def validate_tables(families, family_rows, matrix_rows, node_rows, names):
    require(len(family_rows) == len(families), 'Family-table length mismatch')
    require({r['family'] for r in family_rows} == set(families), 'Family-table IDs mismatch')
    for row in family_rows:
        c = families[row['family']]
        require(integer(row['n_genomes']) == len(c), 'Family occupancy mismatch')
        require(integer(row['n_genes']) == sum(c.values()), 'Family copy count mismatch')
        require(row['class'] == category(len(c)), 'Family category mismatch')
    require(len(matrix_rows) == len(families), 'Matrix row count mismatch')
    require({r['family'] for r in matrix_rows} == set(families), 'Matrix IDs mismatch')
    for row in matrix_rows:
        require(set(row) == {'family', *names.values()}, 'Matrix genome columns mismatch')
        for genome, value in row.items():
            if genome != 'family':
                require(integer(value) == families[row['family']][genome], 'Matrix copy count mismatch')
    require(len(node_rows) == 36 and len({r['node'] for r in node_rows}) == 36, 'Expected 36 unique nodes')
    sets, rebuilt = [], []
    family_sets = [set(c) for c in families.values()]
    for row in node_rows:
        members = row['members'].split('; ')
        require(len(set(members)) == len(members) and set(members) <= set(names.values()), 'Bad clade members')
        require(integer(row['n_genomes']) == len(members), 'Clade size mismatch')
        s = frozenset(members)
        require(s not in sets, 'Duplicate clade membership')
        sets.append(s)
        sizes = [len(f & s) for f in family_sets]
        result = dict(node=integer(row['node']), n_genomes=len(s), members=members,
                      total_families=sum(n > 0 for n in sizes), singleton=sizes.count(1),
                      perfect_core=sizes.count(len(s)), partial=sum(1 < n < len(s) for n in sizes))
        for key in ('total_families', 'singleton', 'perfect_core', 'partial'):
            require(result[key] == integer(row[key]), f'Clade count mismatch: node {row["node"]}, {key}')
        rebuilt.append(result)
    return sets, rebuilt

def compute(data=ROOT / 'data'):
    data = Path(data)
    pin = verify_hashes(data)
    with (data / 'genome_name_map.tsv').open(encoding='utf-8', newline='') as fh:
        mapping = list(csv.reader(fh, delimiter='\t'))
    require(all(len(r) == 2 for r in mapping), 'Malformed genome map')
    require(len({r[0] for r in mapping}) == len(mapping), 'Duplicate genome-map key')
    names = dict(mapping)
    obj = read_json(gzip.decompress((data / 'pangenome.json.gz').read_bytes()).decode('utf-8'))
    families, gene_count = parse_families(obj, names)
    require(len(families) == 21853 and gene_count == 252945, 'Unexpected raw pangenome totals')
    sets, clades = validate_tables(families, rows(data / 'families.tsv'),
                                  rows(data / 'presence_absence.tsv'), rows(data / 'node_table.tsv'), names)
    tree = read_json((data / 'tree_40tips_historical.json').read_text())
    tips, old_clades = tree_clades(tree['tree'])
    new_tips, new_clades = tree_clades((data / 'tree_40tips_corrected.newick').read_text())
    require(tips == new_tips and set(old_clades) == set(new_clades), 'Corrected tree changes clade membership')
    tipnames = {key: names[value['g'][0]] for key, value in tree['ws_refs'].items()}
    require(set(tipnames) == tips and len(set(tipnames.values())) == 40, 'Tree mapping is not bijective')
    require({frozenset(tipnames[t] for t in c) for c in old_clades} == set(sets), 'Node membership differs from tree')
    raw_report = re.split(r'NODE\s+', (data / 'pan_phylo_table.tsv').read_text())[1:]
    require(len(raw_report) == 36, 'Raw report node count mismatch')
    report_ids = set()
    for block in raw_report:
        node_id = int(block.split()[0]); require(node_id not in report_ids, 'Duplicate report node')
        report_ids.add(node_id)
        row = next(r for r in clades if r['node'] == node_id)
        for label, key in [('TOTAL', 'total_families'), ('singleton', 'singleton'), ('partial', 'partial'), ('perfect core', 'perfect_core')]:
            found = re.search(re.escape(label) + r'\s+(\d+)', block)
            require(found is not None and int(found[1]) == row[key], 'Raw report counts differ')
    means = []
    archived_means = rows(data / 'accumulation_curve.tsv')
    for n in sorted({r['n_genomes'] for r in clades}):
        group = [r for r in clades if r['n_genomes'] == n]
        cores = [r['perfect_core'] for r in group]
        result = dict(n_genomes=n, n_clades=len(group), core_mean=mean(cores), core_min=min(cores), core_max=max(cores), pan_mean=mean(r['total_families'] for r in group))
        old = [r for r in archived_means if integer(r['n_genomes']) == n]
        require(len(old) == 1, 'Missing/duplicate accumulation size')
        for key in ('n_clades', 'core_mean', 'core_min', 'core_max', 'pan_mean'):
            require(integer(old[0][key]) == round(result[key]), 'Rounded historical mean mismatch')
        means.append(result)
    require(len(archived_means) == len(means) == 15, 'Unexpected aggregate row count')
    hist = Counter(len(c) for c in families.values())
    cats = Counter(category(len(c)) for c in families.values())
    full = next(r for r in clades if r['n_genomes'] == 40)
    previous = next(r for r in clades if r['n_genomes'] == 39)
    return dict(families=len(families), genes=gene_count, genomes=40,
                occupancy={str(n): hist[n] for n in range(1, 41)}, categories=dict(cats),
                clades=clades, means=means, gain=full['total_families']-previous['total_families'],
                core_loss=previous['perfect_core']-full['perfect_core'],
                added_genome=sorted(set(full['members'])-set(previous['members'])), input_sha256=pin)

def check_destination(path):
    path = Path(path).expanduser().resolve()
    require(not path.exists(), f'Output already exists: {path}')
    for folder in ('data', 'reference', 'scripts', 'tests', 'docs'):
        protected = ROOT / folder
        require(not path.is_relative_to(protected), f'Output inside protected module folder: {folder}')
    return path

def write_tables(result, out):
    (out / 'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for name, values in [('clades.tsv', result['clades']), ('means.tsv', result['means']),
                         ('occupancy.tsv', [dict(n_genomes=n, n_families=result['occupancy'][str(n)], category=category(n)) for n in range(1,41)])]:
        with (out / name).open('w', encoding='utf-8', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=list(values[0]), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            for row in values:
                writer.writerow({k: '; '.join(v) if isinstance(v, list) else v for k,v in row.items()})
