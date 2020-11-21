"""
"""
import numpy as np

try:  # pysal 1.14
    from pysal import weights as pysal_weights
    from pysal.esda import moran
except ImportError:
    try:  # pysal 2.0
        from pysal.lib import weights as pysal_weights
        from pysal.explore.esda import moran
    except ImportError:
        pass


def get_center(matrix, size, span=1):
    """
    returns the average interactionsin the center of the matrix
    """
    k = size // 2 + 1
    kmspan = k - span
    kpspan = k + span + 1
    return sum(matrix[i][j] for i in range(kmspan, kpspan)
               for j in range(kmspan, kpspan)) / (span * 2 + 1)**2


def matrix_to_decay(matrix, size, metric='loop'):
    mid = size // 2
    # loop
    between = [(di + abs(mid - j), yval) 
                   for di, line in ((abs(mid - i), line)
                                     for i, line in enumerate(matrix[:mid]))
                   for j, yval in enumerate(line[:mid])]
    # not loop: upper half
    outside = [(di + abs(mid - j), yval)
                   for di, line in ((abs(mid - i), line)
                                     for i, line in enumerate(matrix))
                   for j, yval in enumerate(line[mid:])]
    # not loop: lower right corner
    outside += [(di + abs(mid - j), yval)
                    for di, line in ((abs(mid - i), line)
                                     for i, line in enumerate(matrix[mid:]))
                    for j, yval in enumerate(line[:mid])]

    if metric == 'normal':
        xvals, yvals = zip(*(between + outside))
        xvals = np.array(xvals)  # distance to center, sum of distance in X and in Y
        yvals = np.array(yvals)  # interactions
        x = size**0.5 - xvals**0.5
        y = (yvals - np.mean(yvals)) / np.std(yvals)
    elif metric == 'loop':
        xvals, yvals = zip(*between)
        x1 = np.array(xvals)  # distance to center, sum of distance in X and in Y
        yvals = np.array(yvals)  # interactions
        y1 = (yvals - np.mean(yvals)) / np.std(yvals)

        xvals, yvals = zip(*outside)
        x2 = np.array(xvals)  # distance to center, sum of distance in X and in Y
        yvals = np.array(yvals)  # interactions
        y2 = (yvals - np.mean(yvals)) / np.std(yvals)

        x = np.concatenate((x1, x2))
        x = size**0.5 - x**0.5
        y = np.concatenate((y1, y2))
    return x, y



def get_weights(matrix, size, width=2, loop=False):
    """
    Computing the weight of each matrix position in relation to the central cell,
    which is the maximum.
    """
    n = {}
    ws = {}
    hwidth = width / 2.
    ihwidth = 1 + 1. / hwidth
    # Going to negative values; degree of correlation + non-correlation
    if loop:
        calculus = lambda x, y: (ihwidth - x / hwidth) * y
    # Using only positive values; degree of correlation
    else:
        calculus = lambda x, y: 1 / x * y
    for i in range(size):
        for j in range(size):
            b = i + j * size
            n[b] = []
            ws[b] = []
            val = matrix[i][j]
            for k in range(-width, width + 1):
                if i + k >= size or i + k < 0:
                    continue
                for l in range(-width, width+1):
                    if j + l >= size or j + l < 0 or k==l==0:
                        continue
                    c = i + k + (j + l) * size
                    n[b].append(c)
                    ws[b].append(calculus(max(abs(k), abs(l)), val))
    return n, ws


def get_MI(matrix, width=2, loop=False, seed=1):
    """
    Computing the Moran Index for each matrix cell
    -
    Moran Index is a measure of multi-dimensional spatial correlation.
    - A positive I value indicates that the feature is surrounded by features with similar values
      -> part of a cluster.
    - A negative I value indicates that the feature is surrounded by features with dissimilar values
      -> an outlier.
    The Local Moran's index can only be interpreted within the context of the computed Z score or p-value.
    """
    mi_stats= {}
    size = len(matrix)
    mi_stats['size'] = size
    matrixlog2 = np.log2(matrix)

    # Using random_seed to block the Moran Index one, and get consistency with all experiments
    np.random.seed(seed)

    n, ws = get_weights(matrix=matrix, size=size, width=width, loop=loop)

    w = pysal_weights.W(n, ws)
    lm = moran.Moran_Local([[matrixlog2[i][j] for i in range(size)]
                                       for j in range(size)],
                                      w, permutations=9999)

    mi_stats["moranI locals"] = lm.p_sim, lm.q

    gm = moran.Moran([[matrixlog2[i][j] for i in range(size)]
                                 for j in range(size)],
                                w, permutations=9999)

    mi_stats["moranI global"] = gm.I, (gm.VI_rand, gm.seI_rand, gm.z_rand, gm.p_rand)

    return mi_stats
