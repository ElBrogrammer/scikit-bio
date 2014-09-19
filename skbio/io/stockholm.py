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

Examples
--------
Assume we have a basic Stockholm file with the following contents::

    # STOCKHOLM 1.0
    seq1         ACC--G-GGGU
    seq2         TCC--G-GGGA
    #=GC SS_cons (((.....)))
    //

>>> from skbio.sequence import RNA
>>> from skbio.alignment import StockholmAlignment
>>> from StringIO import StringIO
>>> sto_in = StringIO("# STOCKHOLM 1.0\\n"
...                   "seq1     ACC--G-GGGU\\nseq2     TCC--G-GGGA\\n"
...                   "#=GC SS_cons (((.....)))\\n//")
>>> sto_records = StockholmAlignment.from_file(sto_in, RNA)
>>> sto = next(sto_records)
>>> print(sto)
# STOCKHOLM 1.0
seq1          ACC--G-GGGU
seq2          TCC--G-GGGA
#=GC SS_cons  (((.....)))
//
>>> sto.gc
{'SS_cons': '(((.....)))'}

We can also write out information by instantiating the StockholmAlignment
object and then printing it.

>>> from skbio.sequence import RNA
>>> from skbio.alignment import StockholmAlignment
>>> seqs = [RNA("ACC--G-GGGU", id="seq1"),
...     RNA("TCC--G-GGGA", id="seq2")]
>>> gf = {
... "RT": ["TITLE1",  "TITLE2"],
... "RA": ["Auth1;", "Auth2;"],
... "RL": ["J Mol Biol", "Cell"],
... "RM": ["11469857", "12007400"]}
>>> sto = StockholmAlignment(seqs, gf=gf)
>>> print(sto)
# STOCKHOLM 1.0
#=GF RN [1]
#=GF RM 11469857
#=GF RT TITLE1
#=GF RA Auth1;
#=GF RL J Mol Biol
#=GF RN [2]
#=GF RM 12007400
#=GF RT TITLE2
#=GF RA Auth2;
#=GF RL Cell
seq1          ACC--G-GGGU
seq2          TCC--G-GGGA
//

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

from collections import defaultdict, OrderedDict

from skbio.alignment import StockholmAlignment
from skbio.sequence import BiologicalSequence
from skbio.io import (register_reader, register_writer, register_sniffer,
                      StockholmFormatError)


@register_sniffer('stockholm')
def _stockholm_sniffer(fh):
    return False, {}


@register_reader('stockholm')
def _stockholm_to_generator(fh, seq_constructor=BiologicalSequence,
                            strict=False):
    return _read_stockholm_records(fh, seq_constructor, strict)


@register_reader('stockholm', StockholmAlignment)
def _stockholm_to_stockholm_alignment(fh, seq_constructor=BiologicalSequence,
                                      strict=False):
    obj = next(_read_stockholm_records(fh, seq_constructor, strict), None)
    if obj is None:
        raise StockholmFormatError("Stockholm file must contain at least one "
                                   "alignment.")
    return obj


@register_writer('stockholm')
def _generator_to_stockholm(gen, fh):
    _write_header(fh)

    for obj in gen:
        _write_stockholm_record(obj, fh, write_header=False)


@register_writer('stockholm', StockholmAlignment)
def _stockholm_alignment_to_stockholm(obj, fh):
    _write_stockholm_record(obj, fh)


def _read_stockholm_records(fh, seq_constructor, strict):
    line = fh.readline()
    if not line.startswith("# STOCKHOLM 1.0"):
        raise StockholmFormatError("Incorrect header found")
    gs_lines = []
    gf_lines = []
    gr_lines = []
    gc_lines = []
    # OrderedDict used so sequences maintain same order as in file
    seqs = OrderedDict()
    for line in fh:
        line = line.strip()
        if line == "" or line.startswith("# S"):
            # skip blank lines or secondary headers
            continue
        elif line == "//":
            # parse the record since we are at its end
            # build the seuence list for alignment construction
            seqs = [seq_constructor(seq, id=_id) for _id, seq in
                    viewitems(seqs)]
            # get length of sequences in the alignment
            seqlen = len(seqs[0][1])

            # parse information lines
            gf = _parse_gf_info(gf_lines)
            gs = _parse_gs_gr_info(gs_lines)
            gr = _parse_gs_gr_info(gr_lines, strict, seqlen)
            gc = _parse_gc_info(gc_lines, strict, seqlen)

            # yield the actual stockholm object
            yield StockholmAlignment(seqs, gf, gs, gr, gc)

            # reset all storage variables
            gs_lines = []
            gf_lines = []
            gr_lines = []
            gc_lines = []
            seqs = OrderedDict()
        # add the metadata lines to the proper lists
        elif line.startswith("#=GF"):
            gf_lines.append(line)
        elif line.startswith("#=GS"):
            gs_lines.append(line)
        elif line.startswith("#=GR"):
            gr_lines.append(line)
        elif line.startswith("#=GC"):
            gc_lines.append(line)
        else:
            lineinfo = line.split()
            # assume sequence since nothing else in format is left
            # in case of interleaved format, need to do check
            if lineinfo[0] in seqs:
                sequence = seqs[lineinfo[0]]
                seqs[lineinfo[0]] = ''.join([sequence, lineinfo[1]])
            else:
                seqs[lineinfo[0]] = lineinfo[1]


def _parse_gf_info(lines):
    """Takes care of parsing GF lines in stockholm plus special cases"""
    parsed = defaultdict(list)
    # needed for making each multi-line RT and NH one string
    rt = []
    nh = []
    lastline = ""
    for line in lines:
        try:
            init, feature, content = line.split(None, 2)
        except ValueError:
            raise StockholmFormatError("Malformed GF line encountered!"
                                       "\n%s" % line.split(None, 2))
        if init != "#=GF":
            raise StockholmFormatError("Non-GF line encountered!")

        # take care of adding multiline RT to the parsed information
        if lastline == "RT" and feature != "RT":
            # add rt line to the parsed dictionary
            rtline = " ".join(rt)
            rt = []
            parsed["RT"].append(rtline)
        elif feature == "RT":
            rt.append(content)
            lastline = feature
            continue

        # Take care of adding multiline NH to the parsed dictionary
        elif lastline == "NH" and feature != "NH":
            nhline = " ".join(nh)
            nh = []
            parsed["NH"].append(nhline)
        elif feature == "NH":
            nh.append(content)
            lastline = feature
            continue

        # add current feature to the parsed information
        parsed[feature].append(content)
        lastline = feature

    if rt:
        rtline = " ".join(rt)
        parsed["RT"].append(rtline)
    if nh:
        nhline = " ".join(nh)
        parsed["NH"].append(nhline)

    # removing unneccessary lists from parsed. Use .items() for py3 support
    for feature, value in parsed.items():
        # list of multi-line features to join into single string if needed
        if feature in ["CC"]:
            parsed[feature] = ' '.join(value)
        elif len(parsed[feature]) == 1:
            parsed[feature] = value[0]
    return parsed


def _parse_gc_info(lines, strict=False, seqlen=-1):
    """Takes care of parsing GC lines in stockholm format"""
    parsed = {}
    for line in lines:
        try:
            init, feature, content = line.split(None, 2)
        except ValueError:
            raise StockholmFormatError("Malformed GC line encountered!\n%s"
                                       % line.split(None, 2))
        if init != "#=GC":
            raise StockholmFormatError("Non-GC line encountered!")

        # add current feature to the parsed information
        if feature in parsed:
            if strict:
                raise StockholmFormatError("Should not have multiple lines "
                                           "with the same feature: %s" %
                                           feature)
        else:
            parsed[feature] = [content]

    # removing unneccessary lists from parsed. Use .items() for py3 support
    for feature, value in parsed.items():
        parsed[feature] = ''.join(value)
        if strict:
            if len(value) != seqlen:
                raise StockholmFormatError("GC must have exactly one char "
                                           "per position in alignment!")

    return parsed


def _parse_gs_gr_info(lines, strict=False, seqlen=-1):
    """Takes care of parsing GS and GR lines in stockholm format"""
    parsed = {}
    parsetype = ""
    for line in lines:
        try:
            init, label, feature, content = line.split(None, 3)
        except ValueError:
            raise StockholmFormatError("Malformed GS/GR line encountered!"
                                       "\n%s" % line.split(None, 3))
        if parsetype == "":
            parsetype = init
        elif init != parsetype:
                raise StockholmFormatError("Non-GS/GR line encountered!")

        if feature in parsed:
            if label in parsed[feature]:
                # interleaved format, so need list of content
                parsed[feature][label].append(content)
            else:
                parsed[feature][label] = [content]
        else:
            parsed[feature] = {label: [content]}

    # join all the crazy lists created during parsing
    for feature in parsed:
        for label, content in parsed[feature].items():
            parsed[feature][label] = ''.join(content)
            if strict:
                if len(parsed[feature][label]) != seqlen:
                    raise StockholmFormatError("GR must have exactly one "
                                               "char per position in the "
                                               "alignment!")
    return parsed


def _write_stockholm_record(obj, fh, write_header=True):
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
