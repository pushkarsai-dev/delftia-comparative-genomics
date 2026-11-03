"""Validate pinned inputs and recompute all family, clade, and plotted counts."""
import argparse
import json
import sys
from figure2_data import ROOT, check_destination, compute, write_tables

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=str, help='New directory for regenerated tables; omit for read-only verification')
    args = parser.parse_args()
    try:
        out = check_destination(args.output_dir) if args.output_dir else None
        result = compute()
        if out is not None:
            out.mkdir(parents=True, exist_ok=False)
            write_tables(result, out)
        print(json.dumps({k:result[k] for k in ('families','genes','genomes','categories','gain','core_loss')}))
    except (ValueError, OSError, KeyError, TypeError, StopIteration) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 2
    return 0

if __name__ == '__main__':
    sys.exit(main())
