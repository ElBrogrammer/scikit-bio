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

from skbio.alignment import StockholmAlignment
from skbio.io import register_reader, register_writer, register_sniffer


@register_sniffer('stockholm')
def _stockholm_sniffer(fh):
    # TODO implement sniffer
    return False, {}


@register_reader('stockholm')
def _stockholm_to_generator(fh):
    pass


@register_reader('stockholm', StockholmAlignment)
def _stockholm_to_stockholm_alignment(fh):
    pass


@register_writer('stockholm')
def _generator_to_stockholm(gen, fh):
    _write_header(fh)

    for obj in gen:
        _write_stockholm(obj, fh, write_header=False)


@register_writer('stockholm', StockholmAlignment)
def _stockholm_alignment_to_stockholm(obj, fh):
    _write_stockholm(obj, fh)


def _write_stockholm(obj, fh, write_header=True):
    # find length of leader info needed to make file pretty. 10 comes from the
    # characters for '#=GR ' and the feature after the sequence label
    spacing = max(map(len, obj.ids())) + 10

    if write_header:
        _write_header(fh)

    _write_gf_markup(fh, obj.gf)
    _write_gs_markup(fh, obj.gs)
    _write_alignment_and_gr_markup(fh, obj, spacing)
    _write_gc_markup(fh, obj.gc, spacing)
    _write_terminator(fh)


def _write_header(fh):
    fh.write('# STOCKHOLM 1.0\n')


def _write_gf_markup(fh, gf):
    # note: everything must be coerced to str in case gf contains ints or
    # floats
    skip_features = {'NH', 'RC', 'RM', 'RN', 'RA', 'RL'}
    for feature, value in gf.items():
        # list of features to skip and parse special later
        if feature in skip_features:
            continue
        # list of features to parse special
        elif feature == 'TN':
            # trees must be in proper order of identifier then tree
            ident = value if isinstance(value, list) else [value]
            tree = gf['NH'] if isinstance(gf['NH'], list) else [gf['NH']]
            for ident, tree in zip(gf['TN'], gf['NH']):
                fh.write(' '.join(['#=GF', 'TN', str(ident)]))
                fh.write('\n')
                fh.write(' '.join(['#=GF', 'NH', str(tree)]))
                fh.write('\n')
        elif feature == 'RT':
            # make sure each reference block stays together
            # set up lists to zip in case some bits are missing
            # create rn list if needed
            default_none = [0]*len(value)
            rn = gf.get('RN', ['[%i]' % x for x in range(1, len(value)+1)])
            rm = gf.get('RM', default_none)
            rt = gf.get('RT', default_none)
            ra = gf.get('RA', default_none)
            rl = gf.get('RL', default_none)
            rc = gf.get('RC', default_none)
            # order: RN, RM, RT, RA, RL, RC
            for n, m, t, a, l, c in zip(rn, rm, rt, ra, rl, rc):
                fh.write(' '.join(['#=GF', 'RN', n]))
                fh.write('\n')
                if m:
                    fh.write(' '.join(['#=GF', 'RM', str(m)]))
                    fh.write('\n')
                if t:
                    fh.write(' '.join(['#=GF', 'RT', str(t)]))
                    fh.write('\n')
                if a:
                    fh.write(' '.join(['#=GF', 'RA', str(a)]))
                    fh.write('\n')
                if l:
                    fh.write(' '.join(['#=GF', 'RL', str(l)]))
                    fh.write('\n')
                if c:
                    fh.write(' '.join(['#=GF', 'RC', str(c)]))
                    fh.write('\n')
        else:
            # normal addition for everything else
            if not isinstance(value, list):
                value = [value]
            for val in value:
                fh.write(' '.join(['#=GF', feature, str(val)]))
                fh.write('\n')


def _write_gs_markup(fh, gs):
    for feature in gs:
        for seqname in gs[feature]:
            fh.write(' '.join(['#=GS', seqname, feature,
                               str(gs[feature][seqname])]))
            fh.write('\n')


def _write_alignment_and_gr_markup(fh, obj, spacing):
    for label, seq in obj.iteritems():
        spacer = ' ' * (spacing - len(label))
        fh.write(spacer.join([label, str(seq)]))
        fh.write('\n')

        # GR markup for sequence
        for feature in viewkeys(obj.gr):
            value = obj.gr[feature][label]
            leaderinfo = ' '.join(['#=GR', label, feature])
            # TODO possible to get a negative length here?
            spacer = ' ' * (spacing - len(leaderinfo))
            fh.write(spacer.join([leaderinfo, value]))
            fh.write('\n')


def _write_gc_markup(fh, gc, spacing):
    for feature, value in viewitems(gc):
        leaderinfo = ' '.join(['#=GC', feature])
        spacer = ' ' * (spacing - len(leaderinfo))
        fh.write(spacer.join([leaderinfo, str(gc[feature])]))
        fh.write('\n')


def _write_terminator(fh):
    fh.write('//\n')
