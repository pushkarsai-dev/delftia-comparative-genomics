"""Independent aggregate expectations, malformed decisions, and overwrite safety."""
import csv,hashlib,json,subprocess,sys,tempfile,unittest,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'scripts/verify_current_analysis.py'
class CurrentAnalysisChecks(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.module=self.root/'module';self.out=self.root/'out'
        results='reconstruction/full_panel_repeat/results'
        for name in ['rbh_calls.tsv','rbh_calls.json','forward_hits.json','reverse_hits.json','environment.json','query_id_check.json','proteome_manifest.json','presence_matrix.tsv']:
            dest=self.module/results/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/results/name,dest)
        for rel in ['data/region002_gene_conservation.csv','data/region002_rbh_matrix.csv']:
            p=self.module/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes((ROOT/rel).read_bytes())
    def run_check(self):return subprocess.run([sys.executable,str(SCRIPT),'--module',str(self.module),'--output',str(self.out)],capture_output=True,text=True)
    def test_real_aggregate_expectations(self):
        result=self.run_check();self.assertEqual(result.returncode,0,result.stderr)
        report=json.loads((self.out/'pair-to-plot-verification.json').read_text())
        self.assertEqual((report['calls'],report['accepted'],report['universal_genes']),(1800,1774,41))
        self.assertEqual(report['variable_counts'],[25,33,37,39]);self.assertEqual(report['contig_distribution']['1'],23)
        self.assertEqual((self.out/'derived_gene_conservation.csv').read_bytes(),(ROOT/'data/region002_gene_conservation.csv').read_bytes())
    def test_duplicate_pair_rejected_without_output(self):
        p=self.module/'reconstruction/full_panel_repeat/results/rbh_calls.tsv'
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        rows[1]=dict(rows[0])
        with p.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');writer.writeheader();writer.writerows(rows)
        self.assertNotEqual(self.run_check().returncode,0);self.assertFalse(self.out.exists())
    def test_contig_mismatch_rejected_without_output(self):
        p=self.module/'data/region002_rbh_matrix.csv'
        with p.open() as f:rows=list(csv.DictReader(f))
        rows[0]['n_contigs']=str(int(rows[0]['n_contigs'])+1)
        with p.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
        self.assertNotEqual(self.run_check().returncode,0);self.assertFalse(self.out.exists())
    def test_existing_output_preserved(self):
        self.out.mkdir();marker=self.out/'keep.txt';marker.write_text('do not overwrite')
        before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in self.module.rglob('*') if p.is_file()}
        self.assertNotEqual(self.run_check().returncode,0)
        self.assertEqual(marker.read_text(),'do not overwrite')
        self.assertEqual(before,{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in self.module.rglob('*') if p.is_file()})

    def mutate_tsv(self, change):
        p=self.module/'reconstruction/full_panel_repeat/results/rbh_calls.tsv'
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        change(rows)
        with p.open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');writer.writeheader();writer.writerows(rows)
        result=self.run_check()
        self.assertNotEqual(result.returncode,0)
        self.assertFalse(self.out.exists())
        return result

    def test_foreign_query_in_rejected_pair_is_rejected(self):
        def change(rows):
            next(v for v in rows if v['reciprocal']=='False')['query']='NOT_IN_45_QUERY_PANEL'
        self.assertIn('decision grid',self.mutate_tsv(change).stderr)

    def test_missing_pair_is_rejected(self):
        self.assertIn('decision grid',self.mutate_tsv(lambda rows: rows.pop()).stderr)

    def test_tsv_only_field_change_is_rejected(self):
        def change(rows):rows[0]['query_len']=str(int(rows[0]['query_len'])+1)
        self.assertIn('field mismatch',self.mutate_tsv(change).stderr)

    def test_altered_raw_hits_are_rejected_without_output(self):
        p=self.module/'reconstruction/full_panel_repeat/results/forward_hits.json'
        data=json.loads(p.read_text());data.pop(next(iter(data)));p.write_text(json.dumps(data))
        result=self.run_check()
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Raw evidence validation failed',result.stderr)
        self.assertFalse(self.out.exists())

    def test_manifests_required(self):
        (self.module/'reconstruction/full_panel_repeat/results/query_id_check.json').unlink()
        self.assertNotEqual(self.run_check().returncode,0)
        self.assertFalse(self.out.exists())

    def test_consistent_contig_falsification_is_rejected(self):
        results=self.module/'reconstruction/full_panel_repeat/results'
        calls=json.loads((results/'rbh_calls.json').read_text())
        with (results/'rbh_calls.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        groups={}
        for row in calls:
            if row['reciprocal']:groups.setdefault(row['genome'],set()).add(row['contig'])
        genome=next(g for g,contigs in groups.items() if len(contigs)==1)
        # Preserve agreement and the plotted one-contig count; only the identity is false.
        for row in calls:
            if row['genome']==genome and row['reciprocal']:row['contig']='NOT_ACTUAL_CONTIG'
        for row in rows:
            if row['genome']==genome and row['reciprocal']=='True':row['contig']='NOT_ACTUAL_CONTIG'
        (results/'rbh_calls.json').write_text(json.dumps(calls))
        with (results/'rbh_calls.tsv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');writer.writeheader();writer.writerows(rows)
        result=self.run_check();self.assertNotEqual(result.returncode,0)
        self.assertIn('contig disagrees',result.stderr);self.assertFalse(self.out.exists())
