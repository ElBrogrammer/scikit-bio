# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

from collections import OrderedDict
from unittest import TestCase, main

from skbio import DNASequence
from skbio.alignment import StockholmAlignment
from skbio.io.stockholm import (_stockholm_alignment_to_stockholm,
                                _generator_to_stockholm)
from skbio.util import get_data_path


class StockholmTests(TestCase):
    def setUp(self):
        seqs = [
            DNASequence("ACC-G-GGTA", id="seq1"),
            DNASequence("TCC-G-GGCA", id="seq2")
        ]
        gf = OrderedDict([
            ("AC", "RF00360"),
            ("BM", ["cmbuild  -F CM SEED",
                    "cmsearch  -Z 274931 -E 1000000"]),
            ("SQ", "9"),
            ("RT", ["TITLE1",  "TITLE2"]),
            ("RN", ["[1]", "[2]"]),
            ("RA", ["Auth1;", "Auth2;"]),
            ("RL", ["J Mol Biol", "Cell"]),
            ("RM", ["11469857", "12007400"]),
            ('RN', ['[1]', '[2]'])
        ])
        gf_trees = OrderedDict({"NH": ["IMATREE", "IMATREETOO"],
                                "TN": ["Tree2", "Tree1"]})
        gs = {"AC": OrderedDict([("seq1", "111"), ("seq2", "222")])}
        gr = {"SS": OrderedDict([("seq1", "1110101111"),
                                 ("seq2", "0110101110")])}
        gc = {"SS_cons": "(((....)))"}

        self.obj_all_markup = StockholmAlignment(
            seqs, gc=gc, gf=gf, gs=gs, gr=gr)

        self.obj_gc_only = StockholmAlignment(
            seqs, gc=gc, gf=None, gs=None, gr=None)

        self.obj_gf_only = StockholmAlignment(
            seqs, gc=None, gf=gf, gs=None, gr=None)

        self.obj_gf_only_trees = StockholmAlignment(
            seqs, gc=None, gf=gf_trees, gs=None, gr=None)

        self.obj_gs_only = StockholmAlignment(
            seqs, gc=None, gf=None, gs=gs, gr=None)

        self.obj_gr_only = StockholmAlignment(
            seqs, gc=None, gf=None, gs=None, gr=gr)

        self.objs = [
            self.obj_all_markup,
            self.obj_gc_only,
            self.obj_gf_only,
            self.obj_gf_only_trees,
            self.obj_gs_only,
            self.obj_gr_only
        ]

        self.fps = map(
            get_data_path,
            ['stockholm_all_markup', 'stockholm_gc_only', 'stockholm_gf_only',
             'stockholm_gf_only_trees', 'stockholm_gs_only',
             'stockholm_gr_only'])

    def test_stockholm_alignment_to_stockholm(self):
        for obj, fp in zip(self.objs, self.fps):
            fh = StringIO()
            _stockholm_alignment_to_stockholm(obj, fh)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_generator_to_stockholm_empty_generator(self):
        # generator that doesn't yield anything
        def empty_generator():
            raise StopIteration()
            yield

        fh = StringIO()
        _generator_to_stockholm(empty_generator(), fh)
        obs = fh.getvalue()
        fh.close()

        exp = '# STOCKHOLM 1.0\n'
        self.assertEqual(obs, exp)

    def test_generator_to_stockholm_single_item(self):
        for obj, fp in zip(self.objs, self.fps):
            def single_item_generator():
                yield obj

            fh = StringIO()
            _generator_to_stockholm(single_item_generator(), fh)
            obs = fh.getvalue()
            fh.close()

            with open(fp, 'U') as fh:
                exp = fh.read()

            self.assertEqual(obs, exp)

    def test_generator_to_stockholm_multiple_items(self):
        def multiple_item_generator():
            for obj in self.obj_all_markup, self.obj_gr_only:
                yield obj

        fh = StringIO()
        _generator_to_stockholm(multiple_item_generator(), fh)
        obs = fh.getvalue()
        fh.close()

        with open(get_data_path('stockholm_all_markup_and_gr_only'),
                  'U') as fh:
            exp = fh.read()

        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
