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
    # find length of leader info needed to make file pretty
    # 10 comes from the characters for '#=GF ' and the feature after label
    infolen = max(map(len, obj.ids())) + 10

    GF_lines = []
    GS_lines = []
    GC_lines = []
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
                    GF_lines.append(' '.join(["#=GF", "TN", str(ident)]))
                    GF_lines.append(' '.join(["#=GF", "NH", str(tree)]))
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
                    GF_lines.append(' '.join(["#=GF", "RN", n]))
                    if m:
                        GF_lines.append(' '.join(["#=GF", "RM", str(m)]))
                    if t:
                        GF_lines.append(' '.join(["#=GF", "RT", str(t)]))
                    if a:
                        GF_lines.append(' '.join(["#=GF", "RA", str(a)]))
                    if l:
                        GF_lines.append(' '.join(["#=GF", "RL", str(l)]))
                    if c:
                        GF_lines.append(' '.join(["#=GF", "RC", str(c)]))
            else:
                # normal addition for everything else
                if not isinstance(value, list):
                    value = [value]
                for val in value:
                    GF_lines.append(' '.join(["#=GF", feature, str(val)]))

    # add GS information if applicable
    if obj.gs:
        for feature in obj.gs:
            for seqname in obj.gs[feature]:
                GS_lines.append(' '.join(["#=GS", seqname, feature,
                                         str(obj.gs[feature][seqname])]))

    # add GC information if applicable
    if obj.gc:
        for feature, value in viewitems(obj.gc):
            leaderinfo = ' '.join(["#=GC", feature])
            spacer = ' ' * (infolen - len(leaderinfo))
            GC_lines.append(spacer.join([leaderinfo,
                                         str(obj.gc[feature])]))

    sto_lines = ["# STOCKHOLM 1.0"] + GF_lines + GS_lines
    # create seq output along with GR info if applicable
    for label, seq in obj.iteritems():
        spacer = ' ' * (infolen - len(label))
        sto_lines.append(spacer.join([label, str(seq)]))
        # GR info added for sequence
        for feature in viewkeys(obj.gr):
            value = obj.gr[feature][label]
            leaderinfo = ' '.join(['#=GR', label, feature])
            spacer = ' ' * (infolen - len(leaderinfo))
            sto_lines.append(spacer.join([leaderinfo, value]))

    sto_lines.extend(GC_lines)
    # add final slashes to end of file
    sto_lines.append('//\n')

    fh.write('\n'.join(sto_lines))
