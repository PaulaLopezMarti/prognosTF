import os

import matplotlib
matplotlib.use('Agg')

from multiprocessing                 import cpu_count
from argparse                        import ArgumentParser
from random                          import getrandbits

from pytadbit.parsers.hic_bam_parser import printime
from pytadbit.utils.file_handling    import mkdir

from meta_waffle.waffle_io           import write_big_matrix, sort_BAMtsv


def main():
    opts         = get_options()
    inbam        = opts.inbam
    resolution   = opts.resolution
    outdir      = opts.outdir
    tmppath      = opts.tmppath
    biases_file  = opts.biases_file
    dry_run      = opts.dry_run
    
    # a bit of hardcoded parameter never hurts
    metric = 'loop'

    printime('Generating huge matrix')
    nheader = write_big_matrix(inbam, resolution, biases_file, outdir, 
                               nchunks=opts.nchunks, wanted_chrom=opts.chrom,
                               wanted_pos1=opts.pos1, wanted_pos2=opts.pos2,
                               dry_run=dry_run, ncpus=opts.ncpus, 
                               tmpdir=tmppath,
                               clean=not opts.dirty, verbose=opts.verbose,
                               waffle_radii=opts.waffle_radii,
                               metric=metric)

    # rand_hash = "%016x" % getrandbits(64)
    # tmpdir = os.path.join(tmppath, '_tmp_%s' % (rand_hash))
    # mkdir(tmpdir)

    #sort all files for only read once per pair of peaks to extract
    # printime('Sorting huge matrix: {}'.format(outfile))
    # sort_BAMtsv(nheader, outfile, tmpdir)

    # os.system('rm -rf {}'.format(tmpdir))

    printime('Done.')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-bam', '--bam', dest='inbam', required=True, default=False,
                        help='Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--out', dest='outdir', required=True, default=False,
                        help='Output directory to store counts')
    parser.add_argument('-b', '--biases', dest='biases_file', default=None,
                        help='Pickle file with biases')
    parser.add_argument('--tmp', dest='tmppath', required=False, default='/tmp',
                        help='[%(default)s] Path to temporary folder')
    parser.add_argument('--chrom', dest='chrom', required=False, default=None,
                        help='Wanted chromosome')
    parser.add_argument('--pos1', dest='pos1', required=False, default=None,
                        type=int, help='Wanted row')
    parser.add_argument('--pos2', dest='pos2', required=False, default=None,
                        type=int, help='Wanted columns')
    parser.add_argument('--keep_tmp', dest='dirty', default=False, action='store_true',
                        help='Keep temporary files for debugging')
    parser.add_argument('--dry_run', dest='dry_run', default=False, action='store_true',
                        help='print job list and exits')
    parser.add_argument('--verbose', dest='verbose', default=False, action='store_true')
    parser.add_argument('-C', dest='ncpus', default=cpu_count(),
                        type=int, help='Number of CPUs used to read BAM')
    parser.add_argument('--nchunks', dest='nchunks', default=100, metavar='INT',
                        type=int, help='''[%(default)s] chunks in which to cut
                        input bam file (in the BAM parsing step)''')
    parser.add_argument('--chunk_size', type=int, default=1000, metavar='INT',
                        help='''[%(default)s]
                        to avoid overloading memory, scans the genomic matrix 
                        in chunks of the given size (in the waffling step)''')
    parser.add_argument('--waffle_radii', type=int, default=10,  metavar='INT',
                        help='''[%(default)s]
                        number of bins around a given position to extract for
                        the waffle.''')

    opts = parser.parse_args()

    return opts


if __name__ == '__main__':
    exit(main())
