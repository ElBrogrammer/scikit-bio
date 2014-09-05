"""
Stockholm format (:mod:`skbio.io.stockholm`)
============================================

.. currentmodule:: skbio.io.stockholm


TODO fill in high-level description

Format Specification
--------------------

TODO fill in format specification

Format Parameters
-----------------

TODO fill in format parameters

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import viewkeys, viewitems

from skbio.alignment import Alignment, StockholmAlignment
from skbio.io import register_reader, register_writer, register_sniffer


@register_sniffer('stockholm')
def _stockholm_sniffer(fh):
    # TODO implement sniffer
    return False, {}


@register_reader('stockholm')
def _stockholm_to_generator(fh):
    pass


@register_reader('stockholm', Alignment)
def _stockholm_to_alignment(fh):
    pass


@register_reader('stockholm', StockholmAlignment)
def _stockholm_to_stockholm_alignment(fh):
    pass


@register_writer('stockholm')
def _generator_to_stockholm(obj, fh):
    pass


@register_writer('stockholm', Alignment)
def _alignment_to_stockholm(obj, fh):
    pass


@register_writer('stockholm', StockholmAlignment)
def _stockholm_alignment_to_stockholm(obj, fh):
    # write header
    fh.write('# STOCKHOLM 1.0\n')

    # find length of leader info needed to make file pretty
    # 10 comes from the characters for '#=GF ' and the feature after label
    infolen = max(map(len, obj.ids())) + 10

    # NOTE: EVERYTHING MUST BE COERCED TO STR in case int or float passed
    # add GF information if applicable
    if obj.gf:
        skipfeatures = set(("NH", "RC", "RM", "RN", "RA", "RL"))
        for feature, value in obj.gf.items():
            # list of features to skip and parse special later
            if feature in skipfeatures:
                continue
            # list of features to parse special
            elif feature == "TN":
                # trees must be in proper order of identifier then tree
                ident = value if isinstance(value, list) else [value]
                tree = obj.gf["NH"] if isinstance(obj.gf["NH"], list) \
                    else [obj.gf["NH"]]
                for ident, tree in zip(obj.gf["TN"], obj.gf["NH"]):
                    fh.write(' '.join(["#=GF", "TN", str(ident)]))
                    fh.write('\n')
                    fh.write(' '.join(["#=GF", "NH", str(tree)]))
                    fh.write('\n')
            elif feature == "RT":
                # make sure each reference block stays together
                # set up lists to zip in case some bits are missing
                # create rn list if needed
                default_none = [0]*len(value)
                rn = obj.gf.get("RN", ["[%i]" % x for x in
                                range(1, len(value)+1)])
                rm = obj.gf.get("RM", default_none)
                rt = obj.gf.get("RT", default_none)
                ra = obj.gf.get("RA", default_none)
                rl = obj.gf.get("RL", default_none)
                rc = obj.gf.get("RC", default_none)
                # order: RN, RM, RT, RA, RL, RC
                for n, m, t, a, l, c in zip(rn, rm, rt, ra, rl, rc):
                    fh.write(' '.join(["#=GF", "RN", n]))
                    fh.write('\n')
                    if m:
                        fh.write(' '.join(["#=GF", "RM", str(m)]))
                        fh.write('\n')
                    if t:
                        fh.write(' '.join(["#=GF", "RT", str(t)]))
                        fh.write('\n')
                    if a:
                        fh.write(' '.join(["#=GF", "RA", str(a)]))
                        fh.write('\n')
                    if l:
                        fh.write(' '.join(["#=GF", "RL", str(l)]))
                        fh.write('\n')
                    if c:
                        fh.write(' '.join(["#=GF", "RC", str(c)]))
                        fh.write('\n')
            else:
                # normal addition for everything else
                if not isinstance(value, list):
                    value = [value]
                for val in value:
                    fh.write(' '.join(["#=GF", feature, str(val)]))
                    fh.write('\n')

    # add GS information if applicable
    if obj.gs:
        for feature in obj.gs:
            for seqname in obj.gs[feature]:
                fh.write(' '.join(["#=GS", seqname, feature,
                                   str(obj.gs[feature][seqname])]))
                fh.write('\n')

    # create seq output along with GR info if applicable
    for label, seq in obj.iteritems():
        spacer = ' ' * (infolen - len(label))
        fh.write(spacer.join([label, str(seq)]))
        fh.write('\n')
        # GR info added for sequence
        for feature in viewkeys(obj.gr):
            value = obj.gr[feature][label]
            leaderinfo = ' '.join(['#=GR', label, feature])
            spacer = ' ' * (infolen - len(leaderinfo))
            fh.write(spacer.join([leaderinfo, value]))
            fh.write('\n')

    # add GC information if applicable
    if obj.gc:
        for feature, value in viewitems(obj.gc):
            leaderinfo = ' '.join(["#=GC", feature])
            spacer = ' ' * (infolen - len(leaderinfo))
            fh.write(spacer.join([leaderinfo, str(obj.gc[feature])]))
            fh.write('\n')

    # add final slashes to end of file
    fh.write('//\n')
