from itertools import combinations, product
import os, sys, girgs, helper
from time import time


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
    """ Generates the edges of bdf-girg"""
    if thr_const_gen is None:
        thr_const_gen = thr_constant
    MMS = [list(tuple) for tuple in bdf.get_min_max_form()]
    MMS_reduced = [list(tuple) for tuple in optimal_min_max_shortening(bdf)]
    return girgs.generateBDFEdges(weights, positions, MMS, MMS_reduced, bdf.get_depth_vol(), thr_constant, thr_const_gen)


def estimate_threshold_constant(bdf, weights, desired_degree, ignore_intersections=True):
    '''
    Estimates the threshold to achieve a certain average degree.

    Args:
        bdf: The boolean distance function
        weights: the weight the graph is going to be sampled with
        desired_degree: The average degree that the graph should have
        ignore_intersections: should intersection between max-terms be considered (recommended)

    Returns:

    '''
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
    '''
    Return the min-max set of a BDF that is shortened to reduce the length as much as possible.
    If less than 2^17 combination exist, the optimal solution is returned.
    Otherwise the described approximation algorithm is use

    Args:
        bdf: The boolean distance function

    Returns: shortened min-max set of the bdf

    '''
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


def step_weight_adjust(positions, weights, min_max_set, thr_con, depth_vol):
    """
    Given positions, weights and a min-max set of the dimensions 
    the positions for the seperated max-girgs are generated and the weight adjusted to have the
    form used by the girg-generator. 
    
    Args:
        positions: The positions the bdf-girg is generated on
        weights: The original weights
        min_max_set: The min-max set of the simplified bdf. All max-sets must have the same length!
        thr_constant: The threshold constant
        thr_exp: The threshold exponent

    Returns: An array of 3-tuples for each max-set. Each tuple contains the adapted positions (1) 
    the scaled weights (2) and the time that was needed for the computation (3).

    """""
    sub_graphs = []
    scaled_weights = helper.weight_scaling(weights, thr_con, depth_vol)
    for i in min_max_set:
        t = time()
        r_positions = helper.filter_by_index(positions, i)
        sub_graphs.append((r_positions, scaled_weights, time() - t))
    return sub_graphs
    
    
def step_girg_gen(gen_list):
    """
        Generates max-girgs from a list of positions and weights

    Args:
        gen_list: Takes a list of tuples, using just the first two entries. The first has to be the positions,
        the second the weights. 

    Returns: A list of tuples containing the generated edges and the time needed for generation.

    """
    result = []
    for i in gen_list:
        t = time()
        edges = girgs.generateEdges(i[1], i[0], float('inf'))
        result.append((edges, time() - t))
    return result
    
    
def step_girg_assemble(gen_list, position, weights, min_max_set, thr_con, depth_vol):
    """
        Assembles a bdf-girg from the generated edges.

    Args:
        gen_list: A list of tuples with the edges in te first entry. Format (edges, -)
        position: positions the graph has been generated with
        weights: weights the graph has been generated with
        min_max_set: min-max-set of the original bdf
        thr_con: the threshold constant
        thr_exp: the threshold exponent

    Returns: The edges of the bdf-girg, and a list of stats for each max-set [number of edges, number of added edges, time for pros]

    """
    edges = set()
    stats = []
    for i in gen_list:
        t = time()
        before_length = len(edges)
        filtered_edges = set(girgs.checkBDFEdges(weights, position, i[0], min_max_set, depth_vol, thr_con));
        edges |= filtered_edges
        stats.append((len(filtered_edges), len(edges) - before_length, time() - t))
    return edges, stats