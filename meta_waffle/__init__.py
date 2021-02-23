#! /usr/bin/env python

"""
"""

from collections import defaultdict
from heapq import heappush, heappop, heappushpop
from gzip import open as gzip_open
from os.path import getsize

from numpy import array


def parse_peaks(cpeaks1, cpeaks2, resolution, first_is_feature, chrom_sizes,
                badcols, section_pos, windows_span):

    def read_line_feature(line):
        '''
        Get information per peak of a feature +/-
        '''
        c, p1, p2, f = line.split()[:4]
        return c, (int(p1) + int(p2)) // 2 // resolution, f

    def read_line_no_feature(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, (int(p1) + int(p2)) // 2 // resolution, ''

    def read_line_no_feature_but(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, (int(p1) + int(p2)) // 2 // resolution, '{}:{}-{}'.format(c, p1, p2)

    peaks1 = open(cpeaks1, "r")
    peaks2 = open(cpeaks2, "r")

    # find out if bed file contain features, or only coordinates
    line = next(peaks1)
    try:
        read_line_feature(line)
        read_line1 = read_line_feature
    except ValueError:
        if first_is_feature:
            read_line1 = read_line_no_feature_but
        else:
            read_line1 = read_line_no_feature
    peaks1.seek(0)

    line = next(peaks2)
    try:
        read_line_feature(line)
        read_line2 = read_line_feature
    except ValueError:
        read_line2 = read_line_no_feature

    max_chrom = dict((c, chrom_sizes[c] // resolution - windows_span)
                     for c in chrom_sizes)

    peaks1.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate1 = set((c, p, f) for c, p, f in map(read_line1, peaks1)
                          if c in max_chrom and windows_span <= p <= max_chrom[c])

    peaks2.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate2 = set((c, p, f) for c, p, f in map(read_line2, peaks2)
                          if c in max_chrom and windows_span <= p <= max_chrom[c])

    peaks1.seek(0)
    npeaks1 = sum(1 for _ in peaks1)
    peaks2.seek(0)
    npeaks2 = sum(1 for _ in peaks2)

    # sort peaks
    bin_coordinate1 = sorted(bin_coordinate1)
    if cpeaks1 == cpeaks2:
        bin_coordinate2 = bin_coordinate1
    else:
        bin_coordinate2 = sorted(bin_coordinate2)

    peaks1.close()
    if cpeaks1 != cpeaks2:
        peaks2.close()

    submatrices = {}
    coord_conv = {}
    bads = set()
    for c, bs, f in bin_coordinate1 + bin_coordinate2:
        pos = section_pos[c][0] + bs
        beg = pos - windows_span
        end = pos + windows_span + 1
        range_ = [(x, p) for x, p in enumerate(range(beg, end))
                  if p not in badcols]
        if not range_:
            bads.add((c, bs, f))
            continue
        submatrices[beg, end] = range_
        coord_conv[c, bs] = beg, end

    bin_coordinate1 = [k for k in bin_coordinate1 if not k in bads]
    bin_coordinate2 = [k for k in bin_coordinate2 if not k in bads]

    return bin_coordinate1, bin_coordinate2, npeaks1, npeaks2, submatrices, coord_conv


def parse_peak_bins(cpeaks1, cpeaks2, resolution, first_is_feature,
                    chrom_sizes, badcols, section_pos, windows_span):

    def read_line_feature(line):
        '''
        Get information per peak of a feature +/-
        '''
        c, p1, p2, f = line.split()[:4]
        return c, (int(p1) + int(p2)) // 2 // resolution, f

    def read_line_no_feature(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, (int(p1) + int(p2)) // 2 // resolution, ''

    def read_line_no_feature_but(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        pos = (int(p1) + int(p2)) // 2 // resolution
        return c, pos, '{}:{}'.format(c, pos)

    peaks1 = open(cpeaks1, "r")
    peaks2 = open(cpeaks2, "r")

    # find out if bed file contain features, or only coordinates
    line = next(peaks1)
    try:
        read_line_feature(line)
        read_line1 = read_line_feature
    except ValueError:
        if first_is_feature:
            read_line1 = read_line_no_feature_but
        else:
            read_line1 = read_line_no_feature
    peaks1.seek(0)

    line = next(peaks2)
    try:
        read_line_feature(line)
        read_line2 = read_line_feature
    except ValueError:
        read_line2 = read_line_no_feature

    max_chrom = dict((c, chrom_sizes[c] // resolution - windows_span)
                     for c in chrom_sizes)

    peaks1.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate1 = set((c, p, f) for c, p, f in map(read_line1, peaks1)
                          if c in max_chrom and windows_span <= p <= max_chrom[c])

    peaks2.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate2 = set((c, p, f) for c, p, f in map(read_line2, peaks2)
                          if c in max_chrom and windows_span <= p <= max_chrom[c])

    peaks1.seek(0)
    npeaks1 = sum(1 for _ in peaks1)
    peaks2.seek(0)
    npeaks2 = sum(1 for _ in peaks2)

    # sort peaks
    bin_coordinate1 = sorted(bin_coordinate1)
    if cpeaks1 == cpeaks2:
        bin_coordinate2 = bin_coordinate1
    else:
        bin_coordinate2 = sorted(bin_coordinate2)

    peaks1.close()
    if cpeaks1 != cpeaks2:
        peaks2.close()

    coord_conv = {}
    bads = set()
    for c, bs, f in bin_coordinate1 + bin_coordinate2:
        pos = section_pos[c][0] + bs
        if pos in badcols:
            bads.add((c, bs))
            continue
        coord_conv[c, bs] = pos

    bin_coordinate1 = [k for k in bin_coordinate1 if not k[:2] in bads]
    bin_coordinate2 = [k for k in bin_coordinate2 if not k[:2] in bads]

    return bin_coordinate1, bin_coordinate2, npeaks1, npeaks2, coord_conv


def generate_pairs(bin_coordinate1, bin_coordinate2, windows_span,
                   window, coord_conv, both_features):

    wsp = (windows_span * 2) + 1

    # put pairs in intervals
    if window == 'inter':
        def test(a, b): return (a[0] != b[0]
                                and a != b)
    elif window == 'intra':
        def test(a, b): return (a[0] == b[0]
                                and wsp <= abs(b[1] - a[1])
                                and a != b)
    elif window == 'all':
        def test(a, b): return ((a[0] == b[0]
                                 and wsp <= abs(b[1] - a[1])
                                 and a != b) or (a[0] != b[0] and a != b))
    else:
        lower, upper = window

        def test(a, b): return (a[0] == b[0]
                                and wsp <= abs(b[1] - a[1])
                                and a != b
                                and lower < abs(a[1] - b[1]) <= upper)

    if bin_coordinate1 is bin_coordinate2:  # we want only one side
        pairs = ((a, b) for i, a in enumerate(bin_coordinate1, 1)
                 for b in bin_coordinate2[i:]
                 if test(a, b))
    else:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if test(a, b))

    # Sort pairs of coordinates according to genomic position of the
    # smallest of each pair, and store it into a new list
    final_pairs = []
    for (chr1, bs1, f1), (chr2, bs2, f2) in pairs:
        beg1, end1 = coord_conv[chr1, bs1]
        beg2, end2 = coord_conv[chr2, bs2]

        what = f1 + f2
        what_new = ''

        if beg1 > beg2:
            if both_features:
                what_new = "{}:{}-{}:{}".format(chr2, bs2, chr1, bs1)
            final_pairs.append((beg2, end2, beg1, end1, what, what_new))

        else:
            if both_features:
                what_new = "{}:{}-{}:{}".format(chr1, bs1, chr2, bs2)
            final_pairs.append((beg1, end1, beg2, end2, what, what_new))

    return sorted(set(final_pairs))


def generate_pair_bins(bin_coordinate1, bin_coordinate2, windows_span,
                       window, coord_conv, first_is_feature, beg_bin, end_bin):

    wsp = (windows_span * 2) + 1

    # put pairs in intervals
    def test(a, b): return ((beg_bin <= coord_conv[a[:2]] <= end_bin or 
                             beg_bin <= coord_conv[b[:2]] <= end_bin)
                            and wsp <= abs(b[1] - a[1])
                            and a != b)
    # if window == 'inter':
    #     def test(a, b): return (a[0] != b[0]
    #                             and a != b)
    # elif window == 'intra':
    #     def test(a, b): return (a[0] == b[0]
    #                             and wsp <= abs(b[1] - a[1])
    #                             and a != b)
    # elif window == 'all':
    #     def test(a, b): return ((a[0] == b[0]
    #                              and wsp <= abs(b[1] - a[1])
    #                              and a != b) or (a[0] != b[0] and a != b))
    # else:
    #     lower, upper = window

    #     def test(a, b): return (a[0] == b[0]
    #                             and wsp <= abs(b[1] - a[1])
    #                             and a != b
    #                             and lower < abs(a[1] - b[1]) <= upper)

    if bin_coordinate1 is bin_coordinate2:  # we want only one side
        pairs = ((a, b) for i, a in enumerate(bin_coordinate1, 1)
                 for b in bin_coordinate2[i:]
                 if test(a, b))
    else:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if test(a, b))

    # Sort pairs of coordinates according to genomic position of the
    # smallest of each pair, and store it into a new list
    if first_is_feature:
        def format_what(a): return a
    else:
        def format_what(a): return ''
    final_pairs = {}
    for (chr1, bs1, f1), (chr2, bs2, f2) in pairs:
        pos1 = coord_conv[chr1, bs1]
        pos2 = coord_conv[chr2, bs2]

        what = f1 + f2
        what_new = format_what(f1)

        if pos1 > pos2:
            # what_new = format_what(chr2, bs2, chr1, bs1)
            final_pairs[pos2, pos1] = what, what_new

        else:
            # what_new = format_what(chr1, bs1, chr2, bs2)
            final_pairs[pos1, pos2] = what, what_new
    return final_pairs


def submatrix_coordinates(final_pairs, wsp, submatrices, counter, both_features):
    '''
    Input BED file(s) of ChIP peaks and bin into desired resolution of Hi-C
    '''
    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features

    # in buf we store a list of coordinates to be yielded
    # when buf spans for twice the window span we sort it and empty it

    buf = []
    for beg1, end1, beg2, end2, what, what_new in final_pairs:
        if both_features:
            counter[''] += 1
        else:
            counter[what] += 1
        range2 = submatrices[beg2, end2]
        for x, p1 in submatrices[beg1, end1]:
            for y, p2 in range2:
                heappush(buf, ((p1, p2), x, y, what, what_new))
        if len(buf) >= wsp:  # need more: genome size times window height
            break

    for beg1, end1, beg2, end2, what, what_new in final_pairs[sum(counter.values()):]:
        if both_features:
            counter[''] += 1
        else:
            counter[what] += 1
        range2 = submatrices[beg2, end2]
        for x, p1 in submatrices[beg1, end1]:
            for y, p2 in range2:
                yield heappushpop(buf, ((p1, p2), x, y, what, what_new))

    while buf:
        yield heappop(buf)


def _update_pos(fh_genome, p):
    """
    update cursor position, and returns corresponding genomic bins
    """
    fh_genome.seek(p)
    _ = next(fh_genome)
    a, b, _ = next(fh_genome).split(None, 2)
    return (int(a), int(b))


def find_previous_line(fh_genome, wanted_pos, initial_position, first_div=2):
    """
    Place the 'cursor' right before a given pair of genomic bins

    :param fh_genome: file handler with genome interactions positionned after
       comments
    :param wanted_pos: couple of coordinates wanted to find
    :param 2 first_div: the closest from top you expect it to be, the largest
       you should set first_div
    """
    prev_pos = initial_position  # we start at least right after comments
    post_pos = getsize(fh_genome.name)
    temp_pos = prev_pos + (post_pos - prev_pos) // first_div

    pos = _update_pos(fh_genome, temp_pos)

    for _ in range(40):  # hardly more than 30 steps (with > 3 billion lines)
        if pos < wanted_pos:
            prev_pos = temp_pos
            temp_pos = (temp_pos + post_pos) // 2
        elif pos > wanted_pos:
            post_pos = temp_pos
            temp_pos = (prev_pos + temp_pos) // 2
        else:
            break
        pos = _update_pos(fh_genome, temp_pos)
    # rewind a bit in case we are in the matching line
    fh_genome.seek(prev_pos)
    l = next(fh_genome)  # and place the cursor at the beginning of aline
    return prev_pos + len(l)


def readfiles(genomic_file, iter_pairs):
    # create empty meta-waffles
    pos = 0
    with open(genomic_file, 'r') as fh1:
        for line in fh1:
            if not line.startswith('#'):
                break
            pos += len(line)
        fh1.seek(pos)

        try:
            pos2, x, y, group, what_new = next(iter_pairs)
            pos = find_previous_line(fh1, pos2, pos)  # place the cursor
            for line in fh1:
                pos += len(line)
                a, b, raw, nrm = line.split('\t')
                pos1 = a, b = int(a), int(b)
                while pos1 >= pos2:
                    if pos1 == pos2:
                        yield pos1, x, y, float(nrm), group, what_new
                        pos2, x, y, group, what_new = next(iter_pairs)
                        if pos1 != pos2:  # some cells in the peak file point to same genomic cell
                            break
                    else:
                        pos2, x, y, group, what_new = next(iter_pairs)
                else:
                    if (a + 2) < pos2[0]:  # relatively big gap
                        pos = find_previous_line(fh1, pos2, pos, first_div=2)
        except StopIteration:
            pass


def interactions_at_intersection(groups, genomic_mat, iter_pairs, submatrices,
                                 bins, window_size, both_features):
    def write_submatrices(X, Y, x, y, nrm, group, what_new):
        c1, b1 = bins[X]
        c2, b2 = bins[Y]
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            c1, b1, c2, b2, x, y, nrm, group, what_new))

    def write_submatrices_both(x, y, nrm, what_new, window_size):
        index = x + y * window_size  # from 2D matrix coordinates to 1D array
        out.write('{}\t{}\t{}\n'.format(what_new, index, nrm))

    readfiles_iterator = readfiles(genomic_mat, iter_pairs)
    if both_features:
        if bins:
            comp_submatrices_path = "{}.gz".format(submatrices)
            out = gzip_open(comp_submatrices_path, 'wt')
            old = ''
            for (X, Y), x, y, nrm, group, what_new in readfiles_iterator:
                groups['']['sum_nrm'][x, y] += nrm
                groups['']['sqr_nrm'][x, y] += nrm**2
                groups['']['passage'][x, y] += 1
                if what_new == old:
                    write_submatrices_both(
                        x, y, round(nrm, 3), '', window_size)
                else:
                    write_submatrices_both(x, y, round(
                        nrm, 3), what_new, window_size)
                    old = what_new
            out.close()
        else:
            for (X, Y), x, y, nrm, group, _ in readfiles_iterator:
                groups['']['sum_nrm'][x, y] += nrm
                groups['']['sqr_nrm'][x, y] += nrm**2
                groups['']['passage'][x, y] += 1
    else:
        if bins:
            out = open(submatrices, 'w')
            for (X, Y), x, y, nrm, group, what_new in readfiles_iterator:
                groups[group]['sum_nrm'][x, y] += nrm
                groups[group]['sqr_nrm'][x, y] += nrm**2
                groups[group]['passage'][x, y] += 1
                write_submatrices(X, Y, x, y, round(nrm, 3), group, what_new)
            out.close()
        else:
            for (X, Y), x, y, nrm, group, _ in readfiles_iterator:
                groups[group]['sum_nrm'][x, y] += nrm
                groups[group]['sqr_nrm'][x, y] += nrm**2
                groups[group]['passage'][x, y] += 1


def interactions_at_intersection_extended_genomic_matrix(
        groups, pair_peaks, genomic_mat, coord_conv):

    # TO DO
    # get FIRST and LAST genomic coordinate of all groups

    first_position = sorted(pair_peaks.keys())[0]
    last_position = sorted(pair_peaks.keys())[-1]
#######
    current_pos = 0
    with open(genomic_mat, 'r') as fh1:
        for line in fh1:
            if not line.startswith('#'):
                break
            current_pos += len(line)
        fh1.seek(current_pos)
    # #    keep the last position
    #     last_pos = find_previous_line(fh1, last_position, current_pos) + 100000
    #     current_pos = find_previous_line(fh1, first_position, current_pos)
    #    go to wanted START position in the file

        # current_pos = find_previous_line(fh1, start, pos)
        for line in fh1:
            current_pos += len(line)
            # if current_pos > last_pos:
            #     break
    #         COPY BELLOW
            a, b, vals = line.split(None, 2)
            pos = a, b = int(a), int(b)
            try:
                group, _ = pair_peaks[pos]
            except KeyError:
                try:
                    group, _ = pair_peaks[pos[::-1]]
                except KeyError:
                    continue
            corr, pval, center, vals = vals.split()
            corr, pval, center = float(corr), float(pval), float(center)
            vals = array([float(v) for v in vals.split(',')])
            groups[group]['sum_nrm'] += vals
            groups[group]['sqr_nrm'] += vals**2
            groups[group]['counter'] += 1


    # for line in fh_genomic_mat:
    #     a, b, vals = line.split(None, 2)
    #     pos = a, b = int(a), int(b)
    #     try:
    #         group, _ = pair_peaks[pos]
    #     except KeyError:
    #         try:
    #             group, _ = pair_peaks[pos[::-1]]
    #         except KeyError:
    #             continue
    #     corr, pval, center, vals = vals.split()
    #     corr, pval, center = float(corr), float(pval), float(center)
    #     vals = array([float(v) for v in vals.split(',')])
    #     groups[group]['sum_nrm'] += vals
    #     groups[group]['sqr_nrm'] += vals**2
    #     groups[group]['counter'] += 1
    fh1.close()
