#!/usr/bin/env python


from random        import random, seed
from collections   import OrderedDict
try:
    from math      import gcd
except ImportError:
    from fractions import gcd
from subprocess    import Popen
from functools     import reduce
from pickle        import dump

from os.path import join  as os_join
from os.path import split as os_split

import numpy as np

from numpy.random  import negative_binomial

from matplotlib    import pyplot as plt


TEST_PATH = os_split(__file__)[0]

if not TEST_PATH:
    TEST_PATH = '.'

QUICK = False

def load_genome(chroms):
    sections = {}
    total = 0
    for c in chroms:
        sections[c] = total
        total += chroms[c]
    bins = {}
    for c in chroms:
        for i in range(chroms[c]):
            bins[sections[c] + i] = c, i
    d = reduce(lambda a, b: gcd(a, b), chroms.values())
    weighted_chroms = [c for c in chroms for _ in range(int(chroms[c]**1.5 / d))]
    return sections, bins, weighted_chroms


def main():
    # some hardcoded defaults of course....
    seed_num = 7
    if QUICK:
        nrnd = 100
    else:
        nrnd = 1000000

    bin_prob = 0.005

    reso = 10000

    if QUICK:
        chroms = OrderedDict([('1', 50), ('2', 30)])
        npeaks = 8
        # probability that an interaction comes from a loop
        loop_prob = 1
    else:
        chroms = OrderedDict([('1', 500), ('2', 300), ('3', 200)])
        npeaks = 20
        # probability that an interaction comes from a loop
        loop_prob = 0.5
    #############################################


    genome_size = sum(chroms.values())

    sections, bins, weighted_chroms = load_genome(chroms)

    seed(seed_num)
    np.random.seed(seed_num)

    cmprts_pos = {}
    bad_cols = {}
    prob = 0.2
    step = 10
    for c in chroms:
        bad_cols[c] = set()
        if not QUICK:
            for _ in range(chroms[c] // 10):
                bad_cols[c].add(int(random() * chroms[c]))
        cmprts_pos[c] = []
        end = 0
        beg = 0
        while end < chroms[c]:
            if random() < prob:
                cmprts_pos[c].append((beg, end))
                beg = end
            end += step
        cmprts_pos[c].append((beg, end))
    cmprts = {}
    for c in cmprts_pos:
        cmprts[c] = {'A': set(p for i, (beg, end) in enumerate(cmprts_pos[c])
                              for p in range(beg, end) if i %2),
                     'B': set(p for i, (beg, end) in enumerate(cmprts_pos[c])
                              for p in range(beg, end) if not i %2)}

    peaks = set()
    peaks1 = set()
    peaks2 = set()
    for c in range(npeaks):
        bin1 = int(random() * (genome_size - 2))
        if random() < 0.4:
            peaks1.add(bin1)
        else:
            peaks2.add(bin1)
        peaks.add(bin1)

    if not QUICK:
        loops = set()
        for bin1 in peaks:
            for bin2 in peaks:
                if random() < 0.1:
                    continue
                if bin1 in peaks1:
                    range1 = 3
                else:
                    range1 = 2
                if bin2 in peaks1:
                    range2 = 3
                else:
                    range2 = 2
                for i in range(range1):
                    for j in range(range2):
                        loops.add((bin1 + i, bin2 + j))
                        loops.add((bin2 + j, bin1 + i))
    else:
        loops = set()
        for bin1 in peaks1:
            for bin2 in peaks2:
                loops.add((bin1, bin2))
                loops.add((bin2, bin1))

    print('generating SAM')
    Popen('mkdir -p {}/tmp'.format(TEST_PATH), shell=True).communicate()
    Popen('mkdir -p {}/data'.format(TEST_PATH), shell=True).communicate()
    out = open(os_join(TEST_PATH, 'tmp', 'fake.sam'), 'w')
    out.write('@HD\tVN:1.5\tSO:coordinate\n')
    for c in chroms:
        out.write('@SQ\tSN:%s\tLN:%d\n' % (c, chroms[c] * reso - 1))
    matrix = [[0 for _ in range(sum(chroms.values()))]
              for _ in range(sum(chroms.values()))]

    nbs = iter(negative_binomial(1, bin_prob, size=nrnd))
    count = 0
    while count < nrnd:
        c1 = weighted_chroms[int(random() * len(weighted_chroms))]
        pos1 = int(random() * chroms[c1])
        if random() > (float(chroms[c1]) / genome_size)**0.8:
            c2 = weighted_chroms[int(random() * len(weighted_chroms))]
            pos2 = int(random() * chroms[c2])
        else:
            c2 = c1
            pos2 = -1
            while pos2 < 0 or pos2 >= chroms[c2]:
                try:
                    if random() < 0.15:
                        wanted_cmprt = 'A' if pos1 in cmprts[c1]['A'] else 'B'
                        while pos2 not in cmprts[c2][wanted_cmprt]:
                            pos2 = pos1 + (next(nbs) * (-1 if random() > 0.5 else 1))
                    else:
                        pos2 = pos1 + (next(nbs) * (-1 if random() > 0.5 else 1))
                except StopIteration:
                    nbs = iter(negative_binomial(1, bin_prob, size=nrnd))
        if pos1 in bad_cols[c1] or pos2 in bad_cols[c2]:
            if random() < 0.5:
                continue
        bin1 = sections[c1] + pos1
        bin2 = sections[c2] + pos2
        if random() <= loop_prob:
            if (bin1, bin2) not in loops:
                continue
        out.write('SRR.{0}\t1024\t{1}\t{2}\t1\t75P\t{3}\t{4}\t75\t*\t*\n'.format(
            count, c1, int(reso / 2 + pos1 * reso), c2, int(reso / 2 + pos2 * reso)))
        out.write('SRR.{0}\t1024\t{1}\t{2}\t1\t75P\t{3}\t{4}\t75\t*\t*\n'.format(
            count, c2, int(reso / 2 + pos2 * reso), c1, int(reso / 2 + pos1 * reso)))
        matrix[bin1][bin2] += 1
        matrix[bin2][bin1] += 1
        count += 1
    out.close()

    print('generating BAM')
    Popen('samtools sort -@ 8 -O BAM {} > {}'.format(
        os_join(TEST_PATH, 'tmp', 'fake.sam'),
        os_join(TEST_PATH, 'data', 'fake.bam')),
          shell=True).communicate()
    Popen('rm -f {}'.format(os_join(TEST_PATH, 'tmp', 'fake.sam')),
          shell=True).communicate()
    Popen('samtools index -@ 8 {}'.format(os_join(TEST_PATH, 'data', 'fake.bam')),
          shell=True).communicate()

    Popen(('tadbit normalize -w {}/tmp --bam {} -r {} --min_count {} '
           '--normalize_only').format(
               TEST_PATH, os_join(TEST_PATH, 'data', 'fake.bam'), reso,
               0 if QUICK else 100), shell=True).communicate()

    Popen(("mv {0}/tmp/04_normalization/biases_*pickle "
           "{0}/data/biases.pickle").format(TEST_PATH), shell=True).communicate()

    Popen('rm -rf {}/tmp'.format(TEST_PATH), shell=True).communicate()

    if QUICK:
        plt.figure(figsize=(10, 7))
    else:
        plt.figure(figsize=(61, 45))
    plt.imshow(np.log2(matrix), interpolation='None', origin='lower')

    total = 0
    xs = []
    ys = []
    for k in chroms:
        xs.append(total)
        xs.append(total)
        ys.append(total)
        total += chroms[k]
        ys.append(total)

    plt.plot(xs, ys, color='k', ls='--')
    plt.plot(ys, xs, color='k', ls='--')
    plt.plot([0, total], [0, total], color='k', alpha=0.5)

    plt.hlines(list(peaks1), 0, list(peaks1), colors='r')
    plt.vlines(list(peaks1), list(peaks1), len(matrix), colors='r')

    plt.hlines(list(peaks2), 0, list(peaks2), colors='b')
    plt.vlines(list(peaks2), list(peaks2), len(matrix), colors='b')

    for p1 in peaks:
        for p2 in peaks:
            if p1 >= p2:
                continue
            if not ((p1 in peaks1 and p2 in peaks2) or
                    (p1 in peaks2 and p2 in peaks1)):
                continue
            for k in range(9):
                for l in range(9):
                    plt.text(p1 + k - 4, p2 + l - 4,
                             matrix[p1 + k - 4][p2 + l - 4],
                             size=2, va='center', ha='center')


    plt.xlim(0, total)
    plt.ylim(0, total)
    plt.colorbar()
    plt.savefig(os_join(TEST_PATH, 'data', 'matrix.pdf'), format='pdf')

    print('saving matrix as pickle')
    out = open(os_join(TEST_PATH, 'data', 'matrix.pickle'), 'wb')
    dump(matrix, out)
    out.close()

    print('Saving BEDs')
    out = open(os_join(TEST_PATH, 'data', 'peaks_protA.bed'), 'w')
    out.write(''.join('{0}\t{1}\t{2}\n'.format(
        bins[p][0],
        bins[p][1] * reso + int(random() * 1000),
        bins[p][1] * reso + int(random() * 1000)) for p in sorted(peaks1)))
    out.close()

    out = open(os_join(TEST_PATH, 'data', 'peaks_protB.bed'), 'w')
    out.write(''.join('{0}\t{1}\t{2}\n'.format(
        bins[p][0],
        bins[p][1] * reso + int(random() * 1000),
        bins[p][1] * reso + int(random() * 1000)) for p in sorted(peaks2)))
    out.close()

    Popen(('cat {0}/data/peaks_protA.bed {0}/data/peaks_protB.bed | '
           'sort -k1n -k2n > {0}/data/peaks_prot.bed').format(TEST_PATH),
          shell=True).communicate()

    out = open(os_join(TEST_PATH, 'data', 'compartments.bed'), 'w')
    out.write(''.join('{}\t{}\t{}\t{}\n'.format(
        c, p * reso, p * reso + reso,
        (1 if p in cmprts[c]['A'] else -1) * (0.2 + 0.8 * random()))
                      for c in chroms
                      for p in range(chroms[c])))
    out.close()

if __name__ == "__main__":
    exit(main())
