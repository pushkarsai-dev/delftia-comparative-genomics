import unittest,sys,csv,json,copy,tempfile,shutil,io,hashlib
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from figure4_data import ROOT,load_data,normalize_track,arrow_vertices,verify_inputs
from rerun_region_search import read_cds
from plot_figure4 import build
from validate_figure4 import check_artists

class Inputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.genes,cls.tracks,cls.totals=load_data()
    def test_hashes(self):verify_inputs(ROOT/'data')
    def test_45_genes(self):self.assertEqual(len(self.genes),45)
    def test_clusterblast_archive(self):
        entries=json.loads((ROOT/'data/knownclusterblast_all.json').read_text())
        self.assertEqual(len(entries),40)
        for entry in entries.values():
            self.assertEqual(hashlib.sha256(entry['text'].encode()).hexdigest(),entry['sha256'])
    def test_20_tracks(self):self.assertEqual(len(self.tracks),20)
    def test_40_panel(self):self.assertEqual(len(self.totals),40)
    def test_reference_order(self):self.assertEqual([g['order'] for g in self.genes],list(range(1,46)))
    def test_20_kcb(self):self.assertEqual(sum(g['kcb'] for g in self.genes),20)
    def test_41_universal(self):self.assertEqual(sum(g['conservation']==40 for g in self.genes),41)
    def test_variable_counts(self):self.assertEqual(sorted(g['conservation'] for g in self.genes if g['conservation']<40),[25,33,37,39])
    def test_conservation_sum(self):self.assertEqual(sum(g['conservation'] for g in self.genes),1774)
    def test_23_single_contig(self):self.assertEqual(sum(int(r['n_contigs'])==1 for r in self.totals),23)
    def test_added_genome(self):self.assertEqual(len(next(t for t in self.tracks if t['genome']=='GCA_037198095.1')['genes']),44)
    def test_original_retained(self):self.assertEqual(len(list(csv.DictReader(io.StringIO((ROOT/'data/region_gene_order.tsv').read_text()),delimiter='\t'))),1293)
    def test_original_exact_rerun(self):
        d=json.loads((ROOT/'data/historical_search_comparison.json').read_text());self.assertEqual(d['missing'],[]);self.assertEqual(d['added'],[])
    def test_cds_identifier_fallback(self):self.assertEqual(len(read_cds(ROOT/'data/regions/GCA_037198095.1.gbk')),44)
    def test_historical_parser_defect(self):self.assertEqual(len(read_cds(ROOT/'data/regions/GCA_037198095.1.gbk','historical')),8)
    def test_second_fallback(self):self.assertEqual(len(read_cds(ROOT/'data/regions/GCA_029835695.1.gbk')),30)
    def test_matched_lengths(self):
        for t in self.tracks:
            for g in t['genes']:self.assertEqual(g['x1']-g['x0'],g['end']-g['start']+1)
    def test_inclusive_span(self):
        for t in self.tracks:self.assertEqual(t['span_bp'],max(g['end'] for g in t['genes'])-min(g['start'] for g in t['genes'])+1)
    def test_hash_corruption_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'input_checksums.json').write_text('{"test":"'+'0'*64+'"}');(p/'test').write_text('bad')
            with self.assertRaises(ValueError):verify_inputs(p)

class Geometry(unittest.TestCase):
    def test_forward_arrow(self):self.assertEqual(arrow_vertices(0,10,0,1)[2],(10,0))
    def test_reverse_arrow(self):self.assertEqual(arrow_vertices(0,10,0,-1)[2],(0,0))
    def test_head_visible(self):self.assertAlmostEqual(arrow_vertices(0,10,0,1)[2][0]-arrow_vertices(0,10,0,1)[1][0],.9)
    def test_short_arrow(self):self.assertAlmostEqual(arrow_vertices(0,.01,0,1)[1][0],.0055)
    def test_bad_interval(self):
        with self.assertRaises(ValueError):arrow_vertices(10,0,0,1)
    def test_bad_strand(self):
        with self.assertRaises(ValueError):arrow_vertices(0,10,0,0)
    def test_reverse_transform(self):
        rows=[dict(start=1,end=100,strand=-1,ref_gene_order=2),dict(start=151,end=200,strand=1,ref_gene_order=1)]
        t,r,s=normalize_track(rows);self.assertTrue(r);self.assertEqual(s,200);self.assertEqual((t[0]['x0'],t[0]['x1'],t[0]['display_strand']),(100,200,1))
    def test_noncollinear_rejected(self):
        rows=[dict(start=i*100+1,end=i*100+70,strand=1,ref_gene_order=o) for i,o in enumerate([1,3,2])]
        with self.assertRaises(ValueError):normalize_track(rows)

class Render(unittest.TestCase):
    def tearDown(self):
        import matplotlib.pyplot as plt
        plt.close('all')
    def test_artist_validation(self):
        fig,_,artists,bars,axb,scale=build();check_artists(fig,artists,bars,axb,scale)
    def test_corrupted_bar_rejected(self):
        fig,_,artists,bars,axb,scale=build();bars[0].set_height(39)
        with self.assertRaises(AssertionError):check_artists(fig,artists,bars,axb,scale)
    def test_wrong_arrow_rejected(self):
        fig,_,artists,bars,axb,scale=build();p=artists[0][2];xy=p.get_xy();xy[:,0]+=1;p.set_xy(xy)
        with self.assertRaises(AssertionError):check_artists(fig,artists,bars,axb,scale)
    def test_wrong_scale_rejected(self):
        fig,_,artists,bars,axb,scale=build();scale.set_position([.3,.3,.4,.01])
        with self.assertRaises(AssertionError):check_artists(fig,artists,bars,axb,scale)

if __name__=='__main__':unittest.main()
