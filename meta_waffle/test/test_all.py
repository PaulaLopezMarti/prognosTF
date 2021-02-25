import unittest
from os.path               import join as os_join
from os.path               import split as os_split
from os import system
from subprocess            import Popen
from collections           import defaultdict
from random                import getrandbits

try:  # python 3
    from pickle            import _Unpickler as Unpickler
except ImportError:  # python 2
    from pickle            import Unpickler

from pytadbit.utils.file_handling    import mkdir

from meta_waffle.utils     import chromosome_from_bam
from meta_waffle           import parse_peaks, generate_pairs
from meta_waffle           import submatrix_coordinates, interactions_at_intersection
from meta_waffle.waffle_io import write_big_matrix, sort_BAMtsv

TEST_PATH = os_join(os_split(os_split(os_split(__file__)[0])[0])[0], 'test')
RESOLUTION = 10000
WINDOWS_SPAN = 4
COORD_CONV = ''
window_size = (WINDOWS_SPAN * 2) + 1
(SECTION_POS, CHROM_SIZES, BINS, PEAK_COORD1, PEAK_COORD2,
 ITER_PAIRS, PAIR_PEAKS) = Unpickler(open(os_join(TEST_PATH, 'test_data.pickle'), 'rb')).load()

GROUPS = Unpickler(open(os_join(TEST_PATH, 'test_result.pickle'), 'rb')).load()


# Popen('python {}/Simulate_HiC.py'.format(TEST_PATH), shell=True).communicate()


class TestWaffle(unittest.TestCase):
    """
    class
    """

    def test_01_bam_loading(self):
        """
        function
        """
        inbam = os_join(TEST_PATH, 'data', 'fake.bam')
        section_pos, chrom_sizes, bins = chromosome_from_bam(
            inbam, RESOLUTION, get_bins=True)
        self.assertEqual(section_pos,
                         {'1': (0, 500), '2': (500, 800), '3': (800, 1000)})
        self.assertEqual(list(chrom_sizes.items()),
                         [('1', 4999999), ('2', 2999999), ('3', 1999999)])
        self.assertEqual(max(bins.keys()), 999)
        self.assertEqual(bins[817], ('3', 17))

    def test_02_parse_peaks(self):
        """
        function
        """
        peak_files = [os_join(TEST_PATH, 'data', 'peaks_protA.bed'),
                      os_join(TEST_PATH, 'data', 'peaks_protB.bed')]
        in_feature = False
        biases = os_join(TEST_PATH, 'data', 'biases.pickle')
        fh = open(biases, "rb")
        try:
            badcols = Unpickler(fh, encoding='latin1').load()['badcol']
        except TypeError:
            badcols = Unpickler(fh).load()['badcol']
        fh.close()
        peak_coord1, peak_coord2, npeaks1, npeaks2, submatrices, coord_conv = parse_peaks(
            peak_files[0], peak_files[1], RESOLUTION, in_feature, CHROM_SIZES, 
            badcols, SECTION_POS, WINDOWS_SPAN)

        global COORD_CONV
        COORD_CONV = coord_conv
        global SUBMATRICES
        SUBMATRICES = submatrices
        self.assertEqual(peak_coord1, PEAK_COORD1)
        self.assertEqual(peak_coord2, PEAK_COORD2)
        self.assertEqual(npeaks1, 6)
        self.assertEqual(npeaks2, 14)

    def test_03_generate_peaks(self):
        """
        function
        """
        max_dist = float('inf')
        window = 'intra'
        pair_peaks = generate_pairs(PEAK_COORD1, PEAK_COORD2,
                                    WINDOWS_SPAN, window, COORD_CONV, both_features=False)
        self.assertEqual(pair_peaks, PAIR_PEAKS)

    def test_04_submatrix_coordinates(self):
        """
        function
        """
        biases = os_join(TEST_PATH, 'data', 'biases.pickle')
        fh = open(biases, "rb")
        try:
            badcols = Unpickler(fh, encoding='latin1').load()['badcol']
        except TypeError:
            badcols = Unpickler(fh).load()['badcol']
        fh.close()
        counter = defaultdict(int)
        iter_pairs = submatrix_coordinates(PAIR_PEAKS, (WINDOWS_SPAN * 2) + 1,
                                           SUBMATRICES, counter, both_features=False)
        iter_pairs = [v for v in iter_pairs]
        self.assertEqual(sorted(iter_pairs), sorted(ITER_PAIRS))
        self.assertEqual(counter[''], 33)


    def test_05_interactions_at_intersection(self):
        """
        function
        """
        genomic_mat = os_join(TEST_PATH, 'data', 'data_bam_10kb.tsv')
        submatrices = os_join(TEST_PATH, 'tmp.tsv')

        groups = {
            '': {
                'sum_nrm' : defaultdict(float),
                'sqr_nrm' : defaultdict(float),
                'passage' : defaultdict(int)
            }
        }

        interactions_at_intersection(
            groups, genomic_mat, (v for v in ITER_PAIRS), submatrices, '',  window_size, both_features=False)
        self.assertEqual(groups['']['passage'], GROUPS['']['passage'])
        self.assertEqual([round(v, 5)for k, v in groups['']['sum_nrm']], 
                         [round(v, 5)for k, v in GROUPS['']['sum_nrm']])
        self.assertEqual([round(v, 5)for k, v in groups['']['sqr_nrm']], 
                         [round(v, 5)for k, v in GROUPS['']['sqr_nrm']])

    def test_06_windows(self):
        """
        test if total intra chromsomal is the same as several windows
        """
        biases = os_join(TEST_PATH, 'data', 'biases.pickle')
        fh = open(biases, "rb")
        try:
            badcols = Unpickler(fh, encoding='latin1').load()['badcol']
        except TypeError:
            badcols = Unpickler(fh).load()['badcol']
        fh.close()
        window = 'intra'
        groups = {}
        windows = [(0, 100), (100, 200), (200, 300), (300, 400)]
        for window in ['intra'] + windows:
            pair_peaks = generate_pairs(PEAK_COORD1, PEAK_COORD2,
                                        WINDOWS_SPAN, window, COORD_CONV, both_features=False)
            counter = defaultdict(int)
            iter_pairs = submatrix_coordinates(pair_peaks, (WINDOWS_SPAN * 1000) + 1, SUBMATRICES, counter, both_features=False)
            genomic_mat = os_join(TEST_PATH, 'data', 'data_bam_10kb.tsv')
            submatrices = os_join(TEST_PATH, 'tmp.tsv')

            groups[window] = {
                '': {
                    'sum_raw' : defaultdict(int),
                    'sqr_raw' : defaultdict(int),
                    'sum_nrm' : defaultdict(float),
                    'sqr_nrm' : defaultdict(float),
                    'passage' : defaultdict(int)
                }
            }

            interactions_at_intersection(groups[window], genomic_mat,
                                         iter_pairs, submatrices, '', window_size, both_features=False)
        self.assertEqual(round(sum(groups['intra']['']['sum_nrm'].values()), 5),
                         round(sum(sum(groups[window]['']['sum_nrm'].values())
                                   for window in windows), 5))
        self.assertEqual(round(sum(groups['intra']['']['sum_nrm'].values()), 5),
                         round(2720.13242866, 5))


    def test_07_big_matrix(self):
        inbam = os_join(TEST_PATH, 'data', 'fake.bam')
        biases = os_join(TEST_PATH, 'data', 'biases3.pickle')
        outfile = os_join(TEST_PATH, 'lele', 'lala.tsv')
        tmppath = os_join(TEST_PATH, 'lele')
        mkdir(tmppath)

        nheader = write_big_matrix(inbam, RESOLUTION, biases, outfile, 
                                   nchunks=100, wanted_chrom=None,
                                   wanted_pos1=None, wanted_pos2=None,
                                   dry_run=False, ncpus=8, 
                                   tmpdir=tmppath,
                                   clean=True, verbose=False,
                                   square_size=100, waffle_radii=WINDOWS_SPAN,
                                   metric='loop')

        rand_hash = "%016x" % getrandbits(64)
        tmpdir = os_join(tmppath, '_tmp_%s' % (rand_hash))
        mkdir(tmpdir)

        #sort all files for only read once per pair of peaks to extract
        sort_BAMtsv(nheader, outfile, tmpdir)

        system('rm -rf {}'.format(tmpdir))
        fh = open(outfile)
        self.assertEqual(187515, sum(1 for l in fh))
        fh.close()
        with open(outfile) as fh:
            for line in fh:
                if line.startswith('525\t723\t'):
                    break
        b, e, r, p, c, vals = line.split()
        self.assertEqual(0.139, float(r))
        self.assertEqual(0.216, float(p))
        self.assertEqual(0.981, float(c))
        self.assertEqual(
            [0.903,0.889,1.401,0.411,0.814,0.417,0.856,0.454,0.8,2.171,0.892,
             4.214,0.433,0,0.402,1.288,0.455,0.869,0.852,0.42,0.919,0,0.842,
             1.579,0.405,1.788,1.706,1.164,1.265,1.328,1.281,1.267,1.249,0,
             0.431,0.428,1.116,1.832,1.698,1.179,0.405,1.996,1.639,0.828,0,
             0.749,0.365,0.383,0.391,1.161,0.795,1.224,0.866,0.786,1.932,1.142,
             1.186,0.732,0.798,0.393,0.421,1.786,0.852,1.366,0.39,0.819,2.621,
             0.741,1.611,0.413,1.371,0.436,1.051,0.345,0,1.165,1.14,0.749,
             1.272,0.45,1.789], [float(v) for v in vals.split(',')])
        with open(outfile) as fh:
            for line in fh:
                if line.startswith('854\t988\t'):
                    break
        b, e, r, p, c, vals = line.split()
        self.assertEqual(0.224, float(r))
        self.assertEqual(0.0448, float(p))
        self.assertEqual(1.394, float(c))
        self.assertEqual(
            [2.123,1.106,0.585,0.572,1.636,1.681,0.517,0.534,0.556,1.521,1.057,
             1.059,1.093,0,0,1.04,1.02,1.062,1.003,2.09,2.093,2.047,1.03,1.058,
             1.028,2.123,0,2.441,1.03,2.062,1.008,1.441,0.521,1.013,0,1.088,
             3.036,1.055,1.069,1.045,0.996,0.512,3.148,2.167,2.256,0.504,2.103,
             0,0,0.993,1.02,1.486,1.621,0.562,0.981,2.044,1.024,0.5,1.447,
             1.983,0.963,1.988,1.639,1.007,1.59,0,0.519,2.473,2.057,1.498,
             0.516,0.537,0.508,1.059,0,0.524,0.499,0.513,0.504,0.521,0],
            [float(v) for v in vals.split(',')])
        system('rm -rf {}'.format(tmppath))

def run():
    unittest.main()


if __name__ == '__main__':
    run()
