"""Full-panel checks, semantic negative controls, and output-safety tests."""
from pathlib import Path
from collections import Counter
import copy
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from figure2_data import (ROOT, category, check_destination, compute, integer,
                         parse_families, read_json, rows, tree_clades, validate_tables, verify_hashes)
from plot_figure2 import create_figure
from validate_figure2 import validate_artists, validate_exported_bars
import matplotlib.pyplot as plt
import pdfplumber

class Figure2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=compute()
        cls.raw=read_json(gzip.decompress((ROOT/'data/pangenome.json.gz').read_bytes()).decode())
        with (ROOT/'data/genome_name_map.tsv').open() as f: cls.names=dict(csv.reader(f,delimiter='\t'))
        cls.families,_=parse_families(cls.raw,cls.names)
        cls.family_rows=rows(ROOT/'data/families.tsv')
        cls.matrix_rows=rows(ROOT/'data/presence_absence.tsv')
        cls.node_rows=rows(ROOT/'data/node_table.tsv')

    def test_all_counts(self):
        self.assertEqual((self.result['families'],self.result['genes'],self.result['genomes']),(21853,252945,40))
        self.assertEqual(self.result['categories'],dict(core=4263,soft_core=307,shell=2576,cloud=14707))

    def test_clades_and_means(self):
        self.assertEqual((len(self.result['clades']),len(self.result['means'])),(36,15))
        self.assertAlmostEqual(self.result['means'][0]['core_mean'],5689.636363636364)
        self.assertEqual(self.result['means'][2]['core_mean'],5387.5)

    def test_endpoints(self):
        self.assertEqual((self.result['gain'],self.result['core_loss']),(312,3))
        self.assertEqual(self.result['added_genome'],['GCA_026410155.1'])
        self.assertEqual((self.result['occupancy']['1'],self.result['occupancy']['40']),(11544,4263))

    def test_category_boundaries(self):
        for n,c in [(1,'cloud'),(5,'cloud'),(6,'shell'),(37,'shell'),(38,'soft_core'),(39,'soft_core'),(40,'core')]:
            with self.subTest(n=n):self.assertEqual(category(n),c)

    def test_invalid_occupancy(self):
        for n in [0,41,-1,True,1.0]:
            with self.subTest(n=n),self.assertRaises(ValueError):category(n)

    def test_invalid_integer(self):
        for n in ['-1','1.0','True','NaN','',True,1]:
            with self.subTest(n=n),self.assertRaises(ValueError):integer(n)

    def test_json_duplicate_key(self):
        with self.assertRaises(ValueError):read_json('{"a":1,"a":2}')

    def test_json_nonfinite(self):
        for n in ['NaN','Infinity','-Infinity']:
            with self.subTest(n=n),self.assertRaises(ValueError):read_json('{"a":'+n+'}')

    def test_duplicate_family(self):
        raw={**self.raw,'orthologs':[self.raw['orthologs'][0],self.raw['orthologs'][0]]}
        with self.assertRaises(ValueError):parse_families(raw,self.names)

    def test_duplicate_gene(self):
        f=copy.deepcopy(self.raw['orthologs'][0]); f['orthologs'].append(f['orthologs'][0])
        with self.assertRaises(ValueError):parse_families({**self.raw,'orthologs':[f]},self.names)

    def test_unknown_genome(self):
        f=copy.deepcopy(self.raw['orthologs'][0]); f['orthologs'][0][2]='unknown'
        with self.assertRaises(ValueError):parse_families({**self.raw,'orthologs':[f]},self.names)

    def test_empty_family(self):
        f={**self.raw['orthologs'][0],'orthologs':[]}
        with self.assertRaises(ValueError):parse_families({**self.raw,'orthologs':[f]},self.names)

    def test_duplicate_map(self):
        names=dict(self.names); keys=list(names); names[keys[0]]=names[keys[1]]
        with self.assertRaises(ValueError):parse_families(self.raw,names)

    def test_wrong_family_category(self):
        rr=[dict(r) for r in self.family_rows]; rr[0]['class']='cloud'
        with self.assertRaises(ValueError):validate_tables(self.families,rr,self.matrix_rows,self.node_rows,self.names)

    def test_wrong_family_count(self):
        rr=[dict(r) for r in self.family_rows]; rr[0]['n_genes']='0'
        with self.assertRaises(ValueError):validate_tables(self.families,rr,self.matrix_rows,self.node_rows,self.names)

    def test_changed_matrix_cell(self):
        rr=[dict(r) for r in self.matrix_rows]; genome=next(iter(self.names.values())); rr[0][genome]='999'
        with self.assertRaises(ValueError):validate_tables(self.families,self.family_rows,rr,self.node_rows,self.names)

    def test_duplicate_matrix_row(self):
        rr=self.matrix_rows[:-1]+[self.matrix_rows[0]]
        with self.assertRaises(ValueError):validate_tables(self.families,self.family_rows,rr,self.node_rows,self.names)

    def test_wrong_clade_count(self):
        rr=[dict(r) for r in self.node_rows]; rr[0]['total_families']='21854'
        with self.assertRaises(ValueError):validate_tables(self.families,self.family_rows,self.matrix_rows,rr,self.names)

    def test_wrong_clade_members(self):
        rr=[dict(r) for r in self.node_rows]; rr[0]['members']+='; unknown'
        with self.assertRaises(ValueError):validate_tables(self.families,self.family_rows,self.matrix_rows,rr,self.names)

    def test_tree_polytomy(self):
        tips,clades=tree_clades('(a:1,(b:1,c:1,d:1):1);')
        self.assertEqual(tips,{'a','b','c','d'});self.assertEqual(len(clades),2)

    def test_duplicate_tree_tip(self):
        with self.assertRaises(ValueError):tree_clades('(a:1,a:1);')

    def test_truncated_tree(self):
        with self.assertRaises(ValueError):tree_clades('(a:1,b:1)')

    def test_plot_artists(self):
        f=create_figure(self.result)
        try:validate_artists(f,self.result)
        finally:plt.close(f)

    def test_wrong_bar_height(self):
        f=create_figure(self.result);f.axes[1].patches[10].set_height(999)
        try:
            with self.assertRaises(ValueError):validate_artists(f,self.result)
        finally:plt.close(f)

    def test_wrong_bar_color(self):
        f=create_figure(self.result);f.axes[1].patches[10].set_facecolor('black')
        try:
            with self.assertRaises(ValueError):validate_artists(f,self.result)
        finally:plt.close(f)

    def test_wrong_mean(self):
        f=create_figure(self.result);y=list(f.axes[0].lines[0].get_ydata());y[-1]+=1;f.axes[0].lines[0].set_ydata(y)
        try:
            with self.assertRaises(ValueError):validate_artists(f,self.result)
        finally:plt.close(f)

    def test_exported_reference_bars_and_negative_controls(self):
        with pdfplumber.open(ROOT/'reference/figure2.pdf') as doc:
            p=doc.pages[0]; bars=sorted([r for r in p.rects if r['height']>1000],key=lambda r:r['x0'])
            axis=next(r for r in p.rects if r['non_stroking_color']==1.0 and r['x0']>p.width/2)
            validate_exported_bars(bars,axis,self.result)
            for field in ['top','x0','non_stroking_color']:
                bad=copy.deepcopy(bars)
                bad[10][field]=(0,0,0) if field=='non_stroking_color' else bad[10][field]+1
                with self.subTest(field=field),self.assertRaises(ValueError):validate_exported_bars(bad,axis,self.result)

    def test_existing_destination(self):
        with tempfile.TemporaryDirectory() as d,self.assertRaises(ValueError):check_destination(d)

    def test_protected_destination(self):
        for d in ['data','reference','scripts','tests','docs']:
            with self.subTest(d=d),self.assertRaises(ValueError):check_destination(ROOT/d/'new-output')

    def test_symlink_destination(self):
        with tempfile.TemporaryDirectory() as d:
            link=Path(d)/'linked';link.symlink_to(ROOT/'data',target_is_directory=True)
            with self.assertRaises(ValueError):check_destination(link/'new-output')

    def test_invalid_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d)/'input_checksums.json').write_text('{}')
            with self.assertRaises(ValueError):verify_hashes(d)

    def test_changed_input_hash(self):
        with tempfile.TemporaryDirectory() as d:
            target=Path(d)/'data';shutil.copytree(ROOT/'data',target)
            with (target/'families.tsv').open('a') as handle:handle.write('\n')
            with self.assertRaisesRegex(ValueError,'Checksum mismatch'):verify_hashes(target)

    def test_caption_and_alt_text(self):
        text=(ROOT/'docs/figure2.md').read_text()
        alt=text.split('## Alt text\n\n',1)[1].split('\n\n',1)[0]
        self.assertLessEqual(len(alt),255)
        for number in ['21,853','4,263','307','2,576','14,707','312']:
            self.assertIn(number,text.split('## Alt text',1)[0])
        self.assertIn('not independent random samples',text)

    def test_cli_refuses_existing_output_without_changes(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'sentinel';p.write_text('unchanged')
            proc=subprocess.run([sys.executable,str(ROOT/'scripts/plot_figure2.py'),'--output-dir',d],capture_output=True,text=True)
            self.assertEqual(proc.returncode,2);self.assertEqual(p.read_text(),'unchanged');self.assertEqual(list(Path(d).iterdir()),[p])

    def test_cli_works_with_optimization(self):
        proc=subprocess.run([sys.executable,'-O',str(ROOT/'scripts/rebuild_counts.py')],capture_output=True,text=True)
        self.assertEqual(proc.returncode,0,proc.stderr)
        self.assertEqual(json.loads(proc.stdout)['families'],21853)

if __name__=='__main__':unittest.main()
