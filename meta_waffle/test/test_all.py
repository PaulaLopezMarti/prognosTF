import unittest
from os.path     import join as os_join
from os.path     import split as os_split
from subprocess  import Popen
from collections import defaultdict

try:  # python 3
    from pickle            import load, _Unpickler as Unpickler
except ImportError:  # python 2
    from pickle            import load, Unpickler

from meta_waffle.utils     import chromosome_from_bam
from meta_waffle.intersect import parse_peaks, generate_pairs
from meta_waffle.intersect import submatrix_coordinates, interactions_at_intersection


TEST_PATH = os_join(os_split(os_split(os_split(__file__)[0])[0])[0], 'test')
RESOLUTION = 10000
WINDOWS_SPAN = 4
window_size = (WINDOWS_SPAN * 2) + 1
(SECTION_POS, CHROM_SIZES, BINS, PEAK_COORD1, PEAK_COORD2,
 ITER_PAIRS, PAIR_PEAKS) = load(open(os_join(TEST_PATH, 'test_data.pickle'), 'rb'))

GROUPS = load(open(os_join(TEST_PATH, 'test_result.pickle'), 'rb'))


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
            peak_files, RESOLUTION, in_feature, CHROM_SIZES, badcols, SECTION_POS, WINDOWS_SPAN, both_features=False)

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
        pair_peaks = generate_pairs(PEAK_COORD1, PEAK_COORD2, RESOLUTION,
                                    WINDOWS_SPAN, max_dist, window, COORD_CONV, both_features=False)
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
        # self.assertEqual(iter_pairs, ITER_PAIRS)
        self.assertEqual(counter[''], 33)


    def test_05_interactions_at_intersection(self):
        """
        function
        """
        genomic_mat = os_join(TEST_PATH, 'data', 'data_bam_10kb.tsv')
        submatrices = os_join(TEST_PATH, 'tmp.tsv')

        groups = {
            '': {
                'sum_raw' : defaultdict(int),
                'sqr_raw' : defaultdict(int),
                'sum_nrm' : defaultdict(float),
                'sqr_nrm' : defaultdict(float),
                'passage' : defaultdict(int)
            }
        }

        interactions_at_intersection(
            groups, genomic_mat, (v for v in ITER_PAIRS), submatrices, '',  window_size, both_features=False)
        self.assertEqual(groups, GROUPS)

    def test_06_windows(self):
        """
        test if total intra chromsomal is the same as several windows
        """
        max_dist = float('inf')
        biases = os_join(TEST_PATH, 'data', 'biases.pickle')
        fh = open(biases, "rb")
        try:
            badcols = Unpickler(fh, encoding='latin1').load()['badcol']
        except TypeError:
            badcols = Unpickler(fh).load()['badcol']
        window = 'intra'
        groups = {}
        windows = [(0, 100), (100, 200), (200, 300), (300, 400)]
        for window in ['intra'] + windows:
            pair_peaks = generate_pairs(PEAK_COORD1, PEAK_COORD2, RESOLUTION,
                                        WINDOWS_SPAN, max_dist, window, COORD_CONV, both_features=False)
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
        self.assertEqual(sum(groups['intra']['']['sum_raw'].values()), 10705)
        self.assertEqual(sum(sum(groups[window]['']['sum_raw'].values())
                             for window in windows), 10705)
        self.assertEqual(round(sum(groups['intra']['']['sum_nrm'].values()), 5),
                         round(sum(sum(groups[window]['']['sum_nrm'].values())
                                   for window in windows), 5))
        self.assertEqual(round(sum(groups['intra']['']['sum_nrm'].values()), 5),
                         round(2720.13242866, 5))


def run():
    unittest.main()


if __name__ == '__main__':
    run()
