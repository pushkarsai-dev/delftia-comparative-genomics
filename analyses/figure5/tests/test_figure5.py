"""Semantic corruption tests bypass checksums to exercise the actual checks."""
import copy, hashlib, json, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from figure5_data import load_data, number, interval
from validate_figure5 import inspect_pdf, sources

class SourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)/'data'
        shutil.copytree(ROOT/'data',self.path)
    def tearDown(self):self.tmp.cleanup()
    def change_json(self,fn):
        path=self.path/'InterProScan_results.json';obj=json.loads(path.read_text());fn(obj);path.write_text(json.dumps(obj))
    def change_tsv(self,fn):
        path=self.path/'InterProScan_results.tsv';rows=[l.split('\t') for l in path.read_text().splitlines()];fn(rows);path.write_text('\n'.join('\t'.join(r) for r in rows)+'\n')
    def reject(self):
        with self.assertRaises((ValueError,KeyError,TypeError)):load_data(self.path,verify_hashes=False)
    def test_real_sources(self):
        p=load_data();self.assertEqual([x['length'] for x in p],[289,329]);self.assertEqual(sum(len(h['fragments']) for x in p for h in x['hits']),9)
    def test_source_read_only(self):
        before={p.name:p.read_bytes() for p in self.path.iterdir()};load_data(self.path);self.assertEqual(before,{p.name:p.read_bytes() for p in self.path.iterdir()})
    def test_checksum(self):
        with (self.path/'InterProScan_results.tsv').open('a') as f:f.write('\n')
        with self.assertRaisesRegex(ValueError,'checksum'):load_data(self.path)
    def test_missing_tsv(self):self.change_tsv(lambda r:r.pop());self.reject()
    def test_duplicate_tsv(self):self.change_tsv(lambda r:r.append(r[0]));self.reject()
    def test_truncated_tsv(self):self.change_tsv(lambda r:r.__setitem__(0,r[0][:5]));self.reject()
    def test_wrong_tsv_evalue(self):self.change_tsv(lambda r:r[0].__setitem__(8,'0.5'));self.reject()
    def test_wrong_tsv_coordinate(self):self.change_tsv(lambda r:r[0].__setitem__(6,'2'));self.reject()
    def test_wrong_tsv_length(self):self.change_tsv(lambda r:r[0].__setitem__(2,'1'));self.reject()
    def test_wrong_tsv_md5(self):self.change_tsv(lambda r:r[0].__setitem__(1,'0'*32));self.reject()
    def test_wrong_tsv_ipr(self):self.change_tsv(lambda r:r[0].__setitem__(11,'IPR000000'));self.reject()
    def test_wrong_tsv_date(self):self.change_tsv(lambda r:r[0].__setitem__(10,'15-08-2026'));self.reject()
    def test_unknown_protein(self):self.change_tsv(lambda r:r[0].__setitem__(0,'unknown'));self.reject()
    def test_missing_json(self):self.change_json(lambda r:r['results'].pop());self.reject()
    def test_duplicate_json(self):self.change_json(lambda r:r['results'].__setitem__(1,r['results'][0]));self.reject()
    def test_wrong_json_sequence(self):self.change_json(lambda r:r['results'][0].__setitem__('sequence','AAA'));self.reject()
    def test_wrong_json_md5(self):self.change_json(lambda r:r['results'][0].__setitem__('md5','0'*32));self.reject()
    def test_missing_signature(self):self.change_json(lambda r:r['results'][0]['matches'].pop());self.reject()
    def test_duplicate_library(self):self.change_json(lambda r:r['results'][0]['matches'].__setitem__(1,r['results'][0]['matches'][0]));self.reject()
    def test_wrong_library_version(self):self.change_json(lambda r:r['results'][0]['matches'][0]['signature']['signatureLibraryRelease'].__setitem__('version','0'));self.reject()
    def test_bad_fragment(self):self.change_json(lambda r:r['results'][0]['matches'][0]['locations'][0]['location-fragments'][0].__setitem__('start',0));self.reject()
    def test_duplicate_fasta(self):
        path=self.path/'BAA554_two_loci.faa';path.write_text(path.read_text()+'\n>HFJBPOEG_02821\nAAAA\n');self.reject()
    def test_invalid_numbers(self):
        for value in (True,False,-1,float('nan'),float('inf'),'1'):
            with self.subTest(value=value),self.assertRaises(ValueError):number(value)
        self.assertEqual(number(0),0)
    def test_invalid_intervals(self):
        for a,b in ((0,2),(3,2),(1,500),(True,2),(1,2.5)):
            with self.subTest(a=a,b=b),self.assertRaises(ValueError):interval(a,b,289)

class ExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.work=Path(cls.tmp.name);cls.out=cls.work/'out'
        cls.before={p.relative_to(ROOT/'data').as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').iterdir()}
        subprocess.run([sys.executable,str(ROOT/'scripts/plot_figure5.py'),'--output-dir',str(cls.out)],cwd=cls.work,check=True,capture_output=True)
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_export(self):
        from validate_figure5 import validate
        self.assertEqual(validate(self.out)['matched_segments'],9)
    def test_no_overwrite(self):
        before={p.name:p.read_bytes() for p in self.out.iterdir()}
        r=subprocess.run([sys.executable,str(ROOT/'scripts/plot_figure5.py'),'--output-dir',str(self.out)],capture_output=True)
        self.assertEqual(r.returncode,2);self.assertEqual(before,{p.name:p.read_bytes() for p in self.out.iterdir()})
    def test_no_source_write(self):
        self.assertEqual(self.before,{p.relative_to(ROOT/'data').as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').iterdir()})
    def test_protected_destination(self):
        r=subprocess.run([sys.executable,str(ROOT/'scripts/plot_figure5.py'),'--output-dir',str(ROOT/'data/test-output')],capture_output=True)
        self.assertEqual(r.returncode,2);self.assertFalse((ROOT/'data/test-output').exists())
    def mutated_pdf_rejected(self,fn,name):
        from plot_figure5 import build_figure
        import matplotlib.pyplot as plt
        p=copy.deepcopy(load_data());fn(p);fig=build_figure(p);path=self.work/(name+'.pdf');fig.savefig(path);plt.close(fig)
        with self.assertRaises(ValueError):inspect_pdf(path,sources(ROOT/'data'))
    def test_filled_gap_rejected(self):
        self.mutated_pdf_rejected(lambda p:p[0]['hits'][1].__setitem__('fragments',[[41,282]]),'filled-gap')
    def test_wrong_segment_rejected(self):
        self.mutated_pdf_rejected(lambda p:p[0]['hits'][0].__setitem__('fragments',[[50,287]]),'wrong-segment')
    def test_wrong_evalue_rejected(self):
        self.mutated_pdf_rejected(lambda p:p[0]['hits'][0].__setitem__('evalue',.001),'wrong-evalue')
    def test_swapped_evalues_rejected(self):
        def swap(p):
            a,b=p[0]['hits'][:2];a['evalue'],b['evalue']=b['evalue'],a['evalue']
        self.mutated_pdf_rejected(swap,'swapped-evalues')
    def test_swapped_accessions_rejected(self):
        def swap(p):
            a,b=p[0]['hits'][:2];a['accession'],b['accession']=b['accession'],a['accession']
        self.mutated_pdf_rejected(swap,'swapped-accessions')
    def test_wrong_motif_rejected(self):
        self.mutated_pdf_rejected(lambda p:p[0]['motif'].update(start=212,end=220),'wrong-motif')
    def test_wrong_length_rejected(self):
        self.mutated_pdf_rejected(lambda p:p[0].__setitem__('length',288),'wrong-length')

if __name__=='__main__':unittest.main()
