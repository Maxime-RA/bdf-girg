import argparse
from time import time
import bdf_girg
import helper
import igraph as ig


def generate_graph(n, deg, ple, bdf, seed):
    positions = bdf_girg.generate_positions(n, len(bdf.get_dimensions()), seed=((seed*3 + 3**9) << 2))
    weights = bdf_girg.generate_weights(n, ple, seed=((seed*2 + 3**9) >> 2))

    thr_con = bdf_girg.estimate_threshold_constant(bdf, weights, deg, config.ignore_intersections)

    if deg < bdf.get_depth_vol() * 4 * bdf.get_length_vol():
        const_high = bdf_girg.estimate_threshold_constant(bdf, weights, bdf.get_depth_vol() * 5 * bdf.get_length_vol(),  config.ignore_intersections)
    else:
        const_high = thr_con

    edges = bdf_girg.gen_bdf_girg(positions, weights, bdf, thr_con, const_high)
    g = ig.Graph(n=n, edges=set(edges))
    del positions, weights, edges
    return g


def generate_grap_time(n, deg, ple, bdf, seed):
    t_pos = time()
    positions = bdf_girg.generate_positions(n, len(bdf.get_dimensions()), seed=((seed * 3 + 3 ** 9) << 2))

    t_wei = time()
    weights = bdf_girg.generate_weights(n, ple, seed=((seed * 2 + 3 ** 9) >> 2))

    t_con = time()
    thr_con = bdf_girg.estimate_threshold_constant(bdf, weights, deg, config.ignore_intersections)

    if deg < bdf.get_depth_vol() * 4 * bdf.get_length_vol():
        const_high = bdf_girg.estimate_threshold_constant(bdf, weights, bdf.get_depth_vol() * 5 * bdf.get_length_vol(),
                                                          config.ignore_intersections)
    else:
        const_high = thr_con

    t_gen = time()
    edges = bdf_girg.gen_bdf_girg(positions, weights, bdf, thr_con, const_high)

    del positions, weights, edges
    return t_wei - t_pos, t_con - t_wei, t_gen - t_con, time() - t_gen


def parse_args():
    # Mandatory arguments
    parser.add_argument('n', type=int, help='Graph size n')
    parser.add_argument('deg', type=float, help='Average degree deg')
    parser.add_argument('ple', type=float, help='Power-law exponent ple')
    parser.add_argument('bdf', type=str, help='Boolean distance function bdf')
    parser.add_argument('seed', type=int, help='Seed for position& weight generation')

    # Optional arguments
    parser.add_argument('-deg_pre', action='store_true', help='Return average degree gotten')
    parser.add_argument('-comp_size', action='store_true', help='Return size of components in descending order')
    parser.add_argument('-deg_dis', action='store_true', help='Returns degree distribution in ascending order')
    parser.add_argument('-dia', action='store_true', help='Return diameter of and effective diameter of graph')
    parser.add_argument('-avg_path_length', action='store_true', help='Return average length of all paths')
    parser.add_argument('-cluster', action='store_true', help='Return average clustering coefficient')
    parser.add_argument('-time', action='store_true', help='Return average clustering coefficient')
    parser.add_argument('-ignore_intersections', action='store_true',
                        help='Ignores the intersections for threshold estimate')

    args = parser.parse_args()

    # Validate mandatory arguments
    if args.n <= 0:
        raise ValueError("Graph must have size of at least 1")
    if not (args.n > args.deg > 1):
        raise ValueError("Average degree must be greater 1 and lower then the number of nodes")
    if not (2.0 < args.ple):
        raise ValueError("Power-law exponent must be greater then 2")

    return args


def get_avg_degree(g):
    return g.ecount() * 2 / g.vcount()


def get_component_size(g):
    return sorted(g.components().sizes(), reverse=True)


def get_degree_dis(g):
    dist = list(g.degree_distribution(bin_width=1).bins())
    result = [0] * int(dist[0][0])
    for (mi, ma, deg) in dist:
        result.append(deg)
        assert len(result) == ma
    return result


def get_diameter(g):
    return g.diameter(directed=False)


def get_average_path_length(g):
    return g.average_path_length()


def get_cluster_coefficient(g):
    return g.transitivity_undirected()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate BDF-GIRG and returning properties')
    config = parse_args()
    # Generate graph
    dst = helper.parse_bdf(config.bdf)
    assert len(dst.get_dimensions()) == max(dst.get_dimensions()) + 1

    if config.time:
        print(','.join(map(str,generate_grap_time(config.n, config.deg, config.ple, helper.parse_bdf(config.bdf), config.seed))))
        exit(0)

    graph = generate_graph(config.n, config.deg, config.ple, helper.parse_bdf(config.bdf), config.seed)

    if config.deg_pre:
        print(get_avg_degree(graph))
    if config.comp_size:
        print(','.join(map(str, get_component_size(graph))))
    if config.deg_dis:
        print(','.join(map(str, get_degree_dis(graph))))
    if config.dia:
        print(get_diameter(graph))
    if config.avg_path_length:
        print(get_average_path_length(graph))
    if config.cluster:
        print(get_cluster_coefficient(graph))
