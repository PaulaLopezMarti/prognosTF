#! /usr/bin/env python
"""
"""

from argparse          import ArgumentParser
try:  # python 3
    from pickle        import _Unpickler as Unpickler, UnpicklingError
except ImportError:  # python 2
    from pickle        import Unpickler, UnpicklingError

from scipy.stats       import spearmanr

from meta_waffle.stats import matrix_to_decay, get_center


def read_pickle(waffle, output, do_loop=True, center_span=1):
    # TODO: deprecated
    out = open(output, 'w')
    try:
        group   = list(waffle.keys())[0]
    except IndexError:
        raise Exception('ERROR: nothing here.')
    size        = waffle[group]['size']
    for group, data in waffle.items():
        counter = data['counter']
        try:
            matrix  = [[data['sum_nrm'][i, j] / counter
                        for i in range(size)[::-1]]
                       for j in range(size)]
            x, y = matrix_to_decay(matrix, size, metric='loop' if do_loop else 'normal')
            spear, pval = spearmanr(x, y)
        except ZeroDivisionError:
            matrix  = [[0 for i in range(size)] for j in range(size)]
            spear = pval = float('nan')
        center = get_center(matrix, size)
        test = pval < 0.05 and spear > 0 and center > 1
        out.write('{}\t{}\t{}\t{}\t{}\n'.format(group, spear, pval, center, test))
    out.close()


def read_text(waffle, output, do_loop=True, center_span=1):
    out = open(output, 'w')
    out.write('# normalize loop: {}\n'.format(
        do_loop))
    out.write('# {}\t{}\t{}\t{}\t{}\n'.format(
        'group name', 'spearman R', 'spearman p-value', 
        'sum interactions in center +-{}'.format(center_span), 
        'test: p-value < 0.05 & center > 1 & R > 0'))
    next(waffle)
    next(waffle)
    next(waffle)
    next(waffle)
    for line in waffle:
        group, resolution, size, counter = line[1:].split()
        size = int(size)
        counter = int(counter)
        if not counter:
            matrix  = [[0 for i in range(size)] for j in range(size)]
            spear = pval = float('nan')
        else:
            vals = (float(n) for n in next(waffle).split())
            matrix = [[next(vals) for _ in range(size)] for _ in range(size)]
            x, y = matrix_to_decay(matrix, size, metric='loop' if do_loop else 'normal')
            spear, pval = spearmanr(x, y)
        center = get_center(matrix, size, span=center_span)
        test = pval < 0.05 and spear > 0 and center > 1
        out.write('{}\t{}\t{}\t{}\t{}\n'.format(group, spear, pval, center, test))
        _ = next(waffle)
    out.close()


def main():
    opts = get_options()

    waffle_file = opts.infile
    output      = opts.outfile
    do_loop     = opts.do_loop
    center_span = opts.center_span

    if opts.is_pickle:
        waffle  = Unpickler(open(waffle_file, 'rb')).load()
        read_pickle(waffle, output, do_loop, center_span)
    else:
        waffle  = open(waffle_file)
        read_text(waffle, output, do_loop, center_span)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='infile', required=True,
                        metavar='PATH', help='path to input file (pickle format)')
    parser.add_argument('-o', dest='outfile', required=True,
                        metavar='PATH', help='''path to output text file with
                        summary statistics by group''')
    parser.add_argument('--loop', dest='do_loop', action='store_true',
                        help='''correct for the formation of
                        tight loops between peaks''')
    parser.add_argument('--center_span', dest='center_span', default=1,
                        type=int,
                        help='''number of bins around center to be considered as center''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')
    parser.add_argument('--pickle', dest='is_pickle', default=False,
                        action='store_true', help='''deprecated''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
