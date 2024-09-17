from itertools import combinations, product
import os, sys, girgs


def generate_weights(n, ple, seed=None):
    """ Generates n weights following a power-law distribution with exponent ple. """
    return girgs.generateWeights(n, ple, seed=seed)


def generate_positions(n, dimension, seed=None):
    """ Generates n random points in the [0,1)^d space. """
    return girgs.generatePositions(n, dimension, seed=seed)


def gen_bdf_edges_trivial(positions, weights, bdf, thr_constant):
    """
    Generates the edges of a bdf-girg using the trivial O(n^2) algorithm.
    Implemented in c++, unpractical for n>10000.
    Main purpose is to test the correct generation on small instances .

    Args:
        positions: Positions in the [0,1)^d space, where d is dimension of the given bdf
        weights: Weights of each position
        bdf: boolean distance function
        thr_constant: threshold exponent controlling the average degree

    Returns: All pairs of vertices (i,j) with bdf distance less or equal to thr_constant * (w[i]*w[j])^(1/D_v(k))

    """
    MMS = [list(tuple) for tuple in bdf.get_min_max_form()]
    return girgs.generateBDFEdgesTrivial(weights, positions, MMS, bdf.get_depth_vol(), thr_constant)


def gen_bdf_girg(positions, weights, bdf, thr_constant, thr_const_gen=None):
    """ Generates the edges of bdf-girg. Basically just makes pre-post- and graph- computation in one step."""
    if thr_const_gen is None:
        thr_const_gen = thr_constant

    MMS = [list(tuple) for tuple in bdf.get_min_max_form()]
    MMS_reduced = [list(tuple) for tuple in optimal_min_max_shortening(bdf)]
    return girgs.generateBDFEdges(weights, positions, MMS, MMS_reduced, bdf.get_depth_vol(), thr_constant, thr_const_gen)


def estimate_threshold_constant(bdf, weights, desired_degree, ignore_intersections=True):
    # Use either simplified volume of all max-girgs or real volume
    if ignore_intersections:
        volume_poly = [0] * (bdf.get_depth_com() + 1)
        for i in bdf.get_min_max_form():
            volume_poly[bdf.get_depth_com() - len(i)] += (2 ** len(i))
    else:
        volume_poly = bdf.get_volume_poly()
        if len(volume_poly) - 1 > 6:
            volume_poly = bdf.get_simplified_poly(1)
    return girgs.estimateThresholdPolynomial(weights, desired_degree, volume_poly,
                                             bdf.get_optimal_bdf(0)[1].get_length_vol(), bdf.get_depth_vol())


def optimal_min_max_shortening(bdf):
    # For every max set, create all possible way of it having size equal to volumetric length
    all_max_set_combi = [set(combinations(t, bdf.get_depth_vol())) for t in list(bdf.get_min_max_form())]

    # For more then 2^17 we use the approximation
    pro = 1
    for i in all_max_set_combi:
        pro *= len(i)
    if pro > 2**17:
        return bdf.get_optimal_bdf(0)[1].get_min_max_form()

    # Try every possible combination
    every_combi = product(*all_max_set_combi)
    min_set = None
    min_length = float('inf')
    for i, combi in enumerate(every_combi):
        s = set(combi)
        if len(s) < min_length:
            min_length = len(s)
            min_set = s
    return min_set
