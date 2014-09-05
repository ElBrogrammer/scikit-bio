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
from skbio.io.stockholm import _stockholm_alignment_to_stockholm


class StockholmAlignmentWriterTests(TestCase):
    def setUp(self):
        self.seqs = [DNASequence("ACC-G-GGTA", id="seq1"),
                     DNASequence("TCC-G-GGCA", id="seq2")]
        self.GF = OrderedDict([
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
        self.GS = {"AC": OrderedDict([("seq1", "111"), ("seq2", "222")])}
        self.GR = {"SS": OrderedDict([("seq1", "1110101111"),
                                      ("seq2", "0110101110")])}
        self.GC = {"SS_cons": "(((....)))"}
        self.st = StockholmAlignment(self.seqs, gc=self.GC, gf=self.GF,
                                     gs=self.GS, gr=self.GR)

    def test_str(self):
        """ Make sure stockholm with all information contained is formatted
        correctly """
        st = StockholmAlignment(self.seqs, gc=self.GC, gf=self.GF, gs=self.GS,
                                gr=self.GR)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GF AC RF00360\n'
               '#=GF BM cmbuild  -F CM SEED\n'
               '#=GF BM cmsearch  -Z 274931 -E 1000000\n'
               '#=GF SQ 9\n'
               '#=GF RN [1]\n'
               '#=GF RM 11469857\n'
               '#=GF RT TITLE1\n'
               '#=GF RA Auth1;\n'
               '#=GF RL J Mol Biol\n'
               '#=GF RN [2]\n'
               '#=GF RM 12007400\n'
               '#=GF RT TITLE2\n'
               '#=GF RA Auth2;\n'
               '#=GF RL Cell\n'
               '#=GS seq1 AC 111\n'
               '#=GS seq2 AC 222\n'
               'seq1          ACC-G-GGTA\n'
               '#=GR seq1 SS  1110101111\n'
               'seq2          TCC-G-GGCA\n'
               '#=GR seq2 SS  0110101110\n'
               '#=GC SS_cons  (((....)))\n//')

        fh = StringIO()
        _stockholm_alignment_to_stockholm(st, fh)
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, exp)

    def test_str_gc(self):
        """ Make sure stockholm with only GC information contained is formatted
        correctly """
        st = StockholmAlignment(self.seqs, gc=self.GC, gf=None, gs=None,
                                gr=None)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n"
               "#=GC SS_cons  (((....)))\n//")

        fh = StringIO()
        _stockholm_alignment_to_stockholm(st, fh)
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, exp)

    def test_str_gf(self):
        """ Make sure stockholm with only GF information contained is formatted
        correctly """
        st = StockholmAlignment(self.seqs, gc=None, gf=self.GF, gs=None,
                                gr=None)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GF AC RF00360\n'
               '#=GF BM cmbuild  -F CM SEED\n'
               '#=GF BM cmsearch  -Z 274931 -E 1000000\n'
               '#=GF SQ 9\n'
               '#=GF RN [1]\n'
               '#=GF RM 11469857\n'
               '#=GF RT TITLE1\n'
               '#=GF RA Auth1;\n'
               '#=GF RL J Mol Biol\n'
               '#=GF RN [2]\n'
               '#=GF RM 12007400\n'
               '#=GF RT TITLE2\n'
               '#=GF RA Auth2;\n'
               '#=GF RL Cell\n'
               'seq1          ACC-G-GGTA\n'
               'seq2          TCC-G-GGCA\n//')

        fh = StringIO()
        _stockholm_alignment_to_stockholm(st, fh)
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, exp)

    def test_str_gs(self):
        """ Make sure stockholm with only GS information contained is formatted
        correctly """
        st = StockholmAlignment(self.seqs, gc=None, gf=None, gs=self.GS,
                                gr=None)
        obs = str(st)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GS seq1 AC 111\n'
               '#=GS seq2 AC 222\n'
               'seq1          ACC-G-GGTA\n'
               'seq2          TCC-G-GGCA\n//')
        self.assertEqual(obs, exp)

    def test_str_gr(self):
        """ Make sure stockholm with only GR information contained is formatted
        correctly """
        st = StockholmAlignment(self.seqs, gc=None, gf=None, gs=None,
                                gr=self.GR)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
               "#=GR seq2 SS  0110101110\n//")

        fh = StringIO()
        _stockholm_alignment_to_stockholm(st, fh)
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, exp)

    def test_str_trees(self):
        """ Make sure stockholm with trees printed correctly"""
        GF = OrderedDict({"NH": ["IMATREE", "IMATREETOO"],
                          "TN": ["Tree2", "Tree1"]})
        st = StockholmAlignment(self.seqs, gc=None, gf=GF, gs=None,
                                gr=None)
        exp = ("# STOCKHOLM 1.0\n#=GF TN Tree2\n#=GF NH IMATREE\n#=GF TN Tree1"
               "\n#=GF NH IMATREETOO\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n//")

        fh = StringIO()
        _stockholm_alignment_to_stockholm(st, fh)
        obs = fh.getvalue()
        fh.close()
        self.assertEqual(obs, exp)


if __name__ == '__main__':
    main()
