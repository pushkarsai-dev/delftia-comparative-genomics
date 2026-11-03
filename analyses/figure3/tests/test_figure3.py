import csv
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from figure3_data import MODULE,FAMILIES,STUDY,EXCEPTION,load_data,reserve_output
from plot_figure3 import build_figure
from validate_figure3 import validate_artists

class Figure3Tests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.data=Path(self.temp.name)/'data'
        shutil.copytree(MODULE/'data',self.data)

    def alter_csv(self, fn):
        path=self.data/'panel40_family_counts.csv'
        with path.open() as h:
            rows=list(csv.DictReader(h))
        fields=list(rows[0])
        fn(rows)
        with path.open('w') as h:
            w=csv.DictWriter(h,fieldnames=fields);w.writeheader();w.writerows(rows)

    def fails(self):
        with self.assertRaises((ValueError,KeyError)):
            load_data(self.data,verify_hashes=False)

    def test_valid(self):
        g,v=load_data();self.assertEqual(np.shape(v),(40,11));self.assertEqual(len(set(g)),40)

    def test_all_original_families_match_fresh(self):
        report=json.loads((self.data/'independent_search.json').read_text())
        with (self.data/'panel40_family_counts.csv').open() as h:
            for r in csv.DictReader(h):
                for k,v in r.items():
                    if k.startswith('PF'):
                        self.assertEqual(int(v),report['counts']['per_genome'][r['genome']][k.split('_')[0]])

    def test_protein_counts(self):
        r=json.loads((self.data/'independent_search.json').read_text())
        self.assertEqual(sum(r['n_per_genome'].values()),239999)

    def test_both_targets_zero(self):
        _,v=load_data();self.assertFalse(np.asarray(v)[:,:2].any())

    def test_one_exception(self):
        g,v=load_data();self.assertEqual([g[i] for i,r in enumerate(v) if r[2]],[EXCEPTION])

    def test_study_first(self):
        g,_=load_data();self.assertEqual(set(g[:4]),set(STUDY))

    def test_duplicate_genome(self):
        self.alter_csv(lambda r:r.__setitem__(1,r[0].copy()));self.fails()

    def test_missing_genome(self):
        self.alter_csv(lambda r:r.pop());self.fails()

    def test_negative_value(self):
        self.alter_csv(lambda r:r[0].__setitem__('PF01977_UbiD','-1'));self.fails()

    def test_empty_value(self):
        self.alter_csv(lambda r:r[0].__setitem__('PF01977_UbiD',''));self.fails()

    def test_fractional_value(self):
        self.alter_csv(lambda r:r[0].__setitem__('PF01977_UbiD','0.5'));self.fails()

    def test_changed_target(self):
        self.alter_csv(lambda r:r[0].__setitem__('PF00282_Pyridoxal_deC','1'));self.fails()

    def test_missing_pf18631_evidence(self):
        (self.data/'pictet_spenglerase_search.txt').write_text('not searched');self.fails()

    def test_missing_fresh_pf18631(self):
        p=self.data/'independent_search.json';r=json.loads(p.read_text())
        del r['counts']['per_genome']['BAA554_hybrid']['PF18631'];p.write_text(json.dumps(r));self.fails()

    def test_altered_fresh_count(self):
        p=self.data/'independent_search.json';r=json.loads(p.read_text())
        r['counts']['per_genome']['BAA554_hybrid']['PF18631']=1;p.write_text(json.dumps(r));self.fails()

    def test_hash_tampering(self):
        with (self.data/'panel40_family_counts.csv').open('a') as h:h.write('\n')
        with self.assertRaises(ValueError):load_data(self.data)

    def test_swapped_counts_preserving_summary(self):
        def swap(rows):
            k='PF02668_TauD';rows[0][k],rows[2][k]=rows[2][k],rows[0][k]
        self.alter_csv(swap);self.fails()

    def test_existing_output_refused(self):
        with self.assertRaises(ValueError):reserve_output(Path(self.temp.name))

    def test_protected_output_refused(self):
        with self.assertRaises(ValueError):reserve_output(MODULE/'data/new-output')

    def test_output_ancestor_refused(self):
        with self.assertRaises(ValueError):reserve_output(MODULE.parent)

    def test_fresh_output(self):
        path=Path(self.temp.name)/'new';self.assertEqual(reserve_output(path),path.resolve())

    def figure(self):
        g,v=load_data();fig,ax=build_figure(g,v);self.addCleanup(plt.close,fig);return g,v,fig,ax

    def test_plot(self):
        g,v,fig,ax=self.figure();validate_artists(fig,ax,g,v)

    def test_wrong_cell_label(self):
        g,v,fig,ax=self.figure();ax.texts[0].set_text('1')
        with self.assertRaises(ValueError):validate_artists(fig,ax,g,v)

    def test_shifted_plot_rows(self):
        g,v,fig,ax=self.figure();ax.images[0].set_data(np.roll(np.asarray(v),1,axis=0))
        with self.assertRaises(ValueError):validate_artists(fig,ax,g,v)

    def test_wrong_genome_label(self):
        g,v,fig,ax=self.figure();labels=[t.get_text() for t in ax.get_yticklabels()];labels[0]='wrong';ax.set_yticklabels(labels)
        with self.assertRaises(ValueError):validate_artists(fig,ax,g,v)

    def test_wrong_colors(self):
        g,v,fig,ax=self.figure();ax.images[0].set_cmap('viridis')
        with self.assertRaises(ValueError):validate_artists(fig,ax,g,v)

    def test_too_small_labels(self):
        g,v,fig,ax=self.figure();ax.texts[0].set_fontsize(4)
        with self.assertRaises(ValueError):validate_artists(fig,ax,g,v)

if __name__=='__main__':
    unittest.main()
