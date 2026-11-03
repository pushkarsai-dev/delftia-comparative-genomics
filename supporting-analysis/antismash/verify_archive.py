#!/usr/bin/env python3
"""Verify the companion antiSMASH archive without extracting or executing it."""
import argparse,hashlib,json,tarfile
from pathlib import Path

def digest_file(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def verify(archive,metadata):
    descriptor=json.loads((metadata/'archive-descriptor.json').read_text())
    entries=json.loads((metadata/'native-output-manifest.json').read_text())
    expected={row['archive_member']:row for row in entries}
    if len(expected)!=len(entries) or len(expected)!=descriptor['member_count']:
        raise ValueError('Manifest member count or uniqueness mismatch')
    if archive.stat().st_size!=descriptor['bytes'] or digest_file(archive)!=descriptor['sha256']:
        raise ValueError('Companion archive byte size or SHA-256 mismatch')
    seen=set();total=0
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            if not member.isfile():raise ValueError('Unexpected non-file archive member')
            if member.name not in expected or member.name in seen:raise ValueError('Unexpected or duplicate archive member')
            row=expected[member.name];seen.add(member.name)
            if member.size!=row['bytes']:raise ValueError('Member size mismatch: '+member.name)
            h=hashlib.sha256();n=0
            with tar.extractfile(member) as handle:
                for chunk in iter(lambda:handle.read(1<<20),b''):h.update(chunk);n+=len(chunk)
            if n!=row['bytes'] or h.hexdigest()!=row['sha256']:raise ValueError('Member checksum mismatch: '+member.name)
            total+=n
    if seen!=set(expected):raise ValueError('Missing archive members')
    return dict(status='PASS',members=len(seen),uncompressed_bytes=total,archive_sha256=descriptor['sha256'])

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(verify(args.archive,Path(__file__).resolve().parent),indent=2))
if __name__=='__main__':main()
