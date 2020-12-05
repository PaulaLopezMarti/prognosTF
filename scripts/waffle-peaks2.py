#! /usr/bin/env python
"""
"""
from meta_waffle import interactions_at_intersection_extended_genomic_matrix
import os
from collections import defaultdict, OrderedDict
from copy        import deepcopy
#from shutil      import copyfileobj

from argparse    import ArgumentParser
try:  # python 3
    from pickle        import dump, HIGHEST_PROTOCOL
except ImportError:  # python 2
    from pickle        import dump, HIGHEST_PROTOCOL

import numpy as np

try:
    from meta_waffle       import parse_peak_bins, generate_pair_bins
    from meta_waffle       import submatrix_coordinates, interactions_at_intersection
    from meta_waffle.utils import printime, mkdir, chromosome_from_header
except ImportError:  # meta-waffle is not installed.. but it's still ok!!!
    from os.path  import join as os_join
    import sys

    sys.path.insert(0, os_join(os.path.split(os.path.split(__file__)[0])[0], 'meta_waffle'))
    from __init__      import parse_peak_bins, generate_pair_bins
    from utils import printime, mkdir, chromosome_from_header


ERROR_INPUT = '''ERROR: file header should be like:
# CHROM	chr1	120000000
# CHROM	chr2	110000000
# CHROM	chr3	90000000
# CHROM	chrX	85000000
# RESOLUTION	50000
# WAFFLE RADII	5
# BADCOLS	1,12,13,54,165,1000
1	7	35	1.2546
2	6	25	2.355
2	8	17	3.0533
2	17	26	1.2519
3	3	11	0.153
...
'''

def main():
    opts = get_options()

    peak_files   = opts.peak_files
    outfile      = opts.outfile
    window       = opts.window
    genomic_mat  = opts.genomic_mat
    in_feature   = opts.first_is_feature
    both_features  = opts.both_are_feature
    submatrix_path  = opts.submatrix_path
    #compress = opts.compress
    silent       = opts.silent

    fh_genomic_mat = open(genomic_mat, 'r')

    chrom_sizes = OrderedDict()
    for line in fh_genomic_mat:
        try:
            rname, chrom, size = line.split('\t')
        except ValueError:
            break
        if not rname.startswith('# CHROM'):
            raise Exception(ERROR_INPUT)
        chrom_sizes[chrom] = int(size)

    if not line.startswith('# RESOLUTION'):
        raise Exception(ERROR_INPUT)

    resolution = int(line.split('\t')[1])

    genome_size = sum(chrom_sizes[c] // resolution + 1 for c in chrom_sizes)

    line = next(fh_genomic_mat)
    if not line.startswith('# WAFFLE RADII'):
        raise Exception(ERROR_INPUT)

    windows_span = int(line.split('\t')[1])

    line = next(fh_genomic_mat)
    if not line.startswith('# BADCOLS'):
        raise Exception(ERROR_INPUT)

    try:
        badcols = set(list(map(int, line.split('\t')[1].split(','))))
    except ValueError:
        badcols = set()

    if window not in  ['inter', 'intra', 'all']:
        window = [int(x) / resolution for x in window.split('-')]
        if window[0] >= window[1]:
            raise Exception('ERROR: beginning of window should be smaller '
                            'than end')

    mkdir(os.path.split(outfile)[0])

    # get chromosome coordinates and converter genomic coordinate to bins
    section_pos, chrom_sizes, bins = chromosome_from_header(
        chrom_sizes, resolution, get_bins=submatrix_path!='')

    # define pairs of peaks
    printime(' - Parsing peaks', silent)
    try:
        peaks1, peaks2 = peak_files
    except ValueError:
        peaks1 = peaks2 = peak_files[0]
    printime(' - Parsing peaks', silent)
    peak_coord1, peak_coord2, npeaks1, npeaks2, coord_conv = parse_peak_bins(
        peaks1, peaks2, resolution, in_feature, chrom_sizes, badcols, 
        section_pos, windows_span)

    # get the groups
    groups = {}
    window_size = (windows_span * 2) + 1

    if both_features:
        groups[''] = {
            'sum_nrm' : np.zeros(window_size**2),
            'sqr_nrm' : np.zeros(window_size**2),
            'passage' : np.zeros(window_size**2)}
        if len(groups) > 1:
            kgroups = list(groups.keys())
            groups = dict(((g1, g2), deepcopy(groups[g1]))
                          for i, g1 in enumerate(kgroups)
                          for g2 in kgroups[i:])
    else:
        for _, _, group in peak_coord1:
            groups[group] = {
                'sum_nrm' : np.zeros(window_size**2),
                'sqr_nrm' : np.zeros(window_size**2),
                'passage' : np.zeros(window_size**2)}

        if not in_feature:
            if len(peak_files) > 1:
                for _, _, group in peak_coord2:
                    groups[group] = {
                        'sum_nrm' : np.zeros(window_size**2),
                        'sqr_nrm' : np.zeros(window_size**2),
                        'passage' : np.zeros(window_size**2)}


    if not silent:
        print((' - Total different (not same bin) and usable (not at chromosome'
               'ends) peaks in {}').format(peak_files[0]))
    printime(('   - {:,} (out of {:,})').format(
        len(peak_coord1), npeaks1), silent)
    if len(peak_files) > 1:
        print((' - Total different (not same bin) and usable (not at chromosome'
               'ends) peaks in {}').format(peak_files[1]))
        printime(('   - {:,} (out of {:,})').format(
            len(peak_coord2), npeaks2), silent)

    printime(' - Generating pairs of coordinates...', silent)
    counter = defaultdict(int)
    pair_peaks = generate_pair_bins(peak_coord1, peak_coord2,
                                    windows_span, window, coord_conv, 
                                    both_features, counter)

    # retrieve interactions at peak pairs using genomic matrix
    # sum them by feature and store them in dictionary
    printime(' - Reading genomic matrix and peaks', silent)
    interactions_at_intersection_extended_genomic_matrix(
        groups, fh_genomic_mat, pair_peaks, both_features)

    printime(' - Submatrices extracted by category:', silent)
    if not silent:
        for group in groups:
            print('    - {:<10} : {:>15}'.format(group if group else 'Total', counter[group]))


    # add the counts of pairs per waffle
    if both_features:
        groups['']['counter']    = counter['']
        groups['']['resolution'] = resolution
        groups['']['size']       = (windows_span * 2) + 1
    else:
        for group in groups:
            groups[group]['counter']    = counter[group]
            groups[group]['resolution'] = resolution
            groups[group]['size']       = (windows_span * 2) + 1

    printime(' - Finished extracting', silent)

    if opts.output_format == 'pickle':
        out = open(outfile, 'wb')
        dump(groups, out, protocol=HIGHEST_PROTOCOL)
        out.close()
    else:
        out = open(outfile, 'w')
        out.write('# Waffle-peak output file\n')
        out.write('# >group name\tresolution\tsub-matrix size\tnumber of submatrices\n')
        out.write('# Sum of normalized interactions\n')
        out.write('# Sum of square normalized interactions\n')
        for group in groups:
            size = groups[group]['size']
            out.write('>{}\t{}\t{}\t{}\n'.format(
                group, groups[group]['resolution'], size, 
                groups[group]['counter']))
            out.write('{}\n'.format(
                '\t'.join(str(round(groups[group]['sum_nrm'][i], 3)) 
                          for i in range(size**2))))
            out.write('{}\n'.format(
                '\t'.join(str(round(groups[group]['sqr_nrm'][i], 3)) 
                          for i in range(size**2))))
        out.close()

    printime('all done!', silent)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('--peaks', dest='peak_files', required=True,
                        nargs="+", metavar='PATH',
                        help='''one or two pairwise peaks files to
                        compute average submatrix (norm and raw). These files
                        should contain at least two columns, chromosome and
                        position, and may contain an extra column of feature. If
                        present, the result will be returned according to the
                        possible combination of this feature''')
    parser.add_argument('--genomic_matrix', dest='genomic_mat', default=True,
                        metavar='GENOMIC_MATRIX', required=True,
                        help='''Path to genomic matrix in 3 columns format
                        (should be sorted with `sort -k1,2n`)''')
    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        metavar='PATH', help='path to output file')
    parser.add_argument('--output_format', default='tsv', choices=['tsv', 'pickle'],
                        metavar='PATH', help='''[%(default)s] path to output file (available
                        formatsare: %(choices)s) in tsv format, parse matrix 
                        like this: 
                        zip([(i, j) for i in range(size) for j in range(size)], line.split())
                        ''')
    parser.add_argument('--all_submatrices', dest='submatrix_path', default='',
                        metavar='PATH', help='''if PATH is provided here, stores
                        all the individual submatrices generated''')
    parser.add_argument('-w', '--window', dest='window', required=False,
                        default='intra', metavar='INT-INT', type=str,
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000"
                        correspond to the window interval, from 1Mb to 2Mb.
                        Use "-w inter" for inter-chromosomal regions, "-w intra" for
                        intra-chromosomal, "-w all" for all combinations
                        (without distance restriction)''')
    parser.add_argument('--first_is_feature', dest='first_is_feature', default=False,
                        action='store_true', help='''When 2 BED files are input,
                        the peaks in the first BED should be also considered as
                        feature. This is to create average sub-matrices for
                        each peak in the first BED file.''')
    parser.add_argument('--both_are_feature', dest='both_are_feature', default=False,
                        action='store_true', help='''When 2 BED files are input,
                        both BED peaks should be considered as feature.
                        Together with --all_submatrices, it creates an average sub-matrices
                        for each pair of peaks and outputs a compress file ''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
