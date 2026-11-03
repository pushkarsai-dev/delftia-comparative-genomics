#!/usr/bin/env python3
"""Compare a full-panel repeat with the deposited reconstruction, read-only.

Require identical selected pairs, decisions, sequences, and matrix. Report,
rather than hide, numerical and lower-ranked alignment differences.
Run validate_evidence.py on both directories separately for raw-hit validity.
"""
import argparse
import collections
import json
import math
from pathlib import Path


def compare(reference, repeat):
    def load(base, name):
        return json.loads((base / name).read_text())

    def index(rows):
        result = {(r['genome'], r['query']): r for r in rows}
        if len(result) != len(rows):
            raise ValueError('Duplicate gene-genome pair')
        return result

    old, new = index(load(reference, 'rbh_calls.json')), index(load(repeat, 'rbh_calls.json'))
    if old.keys() != new.keys() or len(new) != 1800:
        raise ValueError('Expected identical 1,800 gene-genome pairs')
    numerical = {'forward_evalue', 'reverse_evalue', 'forward_score', 'reverse_score'}
    changes = collections.Counter()
    maximum = collections.defaultdict(float)
    relative_e = 0.0
    for pair, current in new.items():
        previous = old[pair]
        if previous.keys() != current.keys():
            raise ValueError(f'Call schema differs: {pair}')
        for field in current:
            a, b = previous[field], current[field]
            if a != b:
                if field not in numerical:
                    raise ValueError(f'Non-numerical call difference: {pair}, {field}')
                if not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in (a, b)):
                    raise ValueError(f'Invalid numerical difference: {pair}, {field}')
                changes[field] += 1
                maximum[field] = max(maximum[field], abs(a-b))
                if field.endswith('evalue') and a > 0:
                    relative_e = max(relative_e, abs(b-a)/a)
    for name in ('query_id_check.json', 'proteome_manifest.json'):
        if load(reference, name) != load(repeat, name):
            raise ValueError(f'Protein/query evidence differs: {name}')
    if (reference/'presence_matrix.tsv').read_bytes() != (repeat/'presence_matrix.tsv').read_bytes():
        raise ValueError('Presence matrix differs')
    raw = {}
    for name in ('forward_hits.json', 'reverse_hits.json'):
        a, b = load(reference, name), load(repeat, name)
        if a.keys() != b.keys():
            raise ValueError(f'Genome set differs in {name}')
        changed = collections.Counter()
        ordering, coordinates, hits, maximum_score = [], [], 0, 0.0
        for genome in a:
            if a[genome].keys() != b[genome].keys():
                raise ValueError(f'Query set differs in {name}: {genome}')
            for query, left in a[genome].items():
                right = b[genome][query]
                lm, rm = ({h['target']:h for h in rows} for rows in (left, right))
                if len(lm) != len(left) or len(rm) != len(right) or lm.keys() != rm.keys():
                    raise ValueError(f'Hit target set differs: {name}, {genome}, {query}')
                if list(lm) != list(rm):
                    ordering.append({'genome':genome, 'query':query})
                for target, u in lm.items():
                    v = rm[target]
                    if u.keys() != v.keys():
                        raise ValueError('Raw-hit schema differs')
                    hits += 1
                    for field in u:
                        if u[field] != v[field]:
                            changed[field] += 1
                            if field in ('score','dom_score'):
                                maximum_score = max(maximum_score, abs(u[field]-v[field]))
                            if field not in {'evalue','score','bias','dom_evalue','dom_score'}:
                                coordinates.append({'genome':genome,'query':query,'target':target,
                                    'field':field,'reference':u[field],'repeat':v[field],
                                    'selected_best_hit':bool(left and target == left[0]['target'])})
        raw[name] = {'hits_compared':hits, 'identical_target_sets':True,
                     'changed_fields':dict(changed),'changed_order_lists':ordering,
                     'coordinate_differences':coordinates,'max_bit_score_difference':maximum_score}
    total = sum(c['reciprocal'] for c in new.values())
    if total != 1774:
        raise ValueError(f'Unexpected reciprocal total: {total}')
    return {'call_count':len(new),'reciprocal_hits':total,'identical_selected_pairs_and_decisions':True,
            'presence_matrix_byte_identical':True,'protein_manifests_equal':True,
            'changed_call_fields':dict(changes),'max_absolute_call_differences':dict(maximum),
            'max_relative_call_evalue_difference':relative_e,'raw_hit_comparisons':raw,
            'scope':'Repeat of the deposited reconstruction, not recovered historical tables.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--repeat', type=Path, required=True)
    args = p.parse_args()
    try:
        print(json.dumps(compare(args.reference, args.repeat), indent=2))
    except (ValueError, KeyError, TypeError, OSError) as error:
        p.exit(2, f'REFUSED: {error}\n')
