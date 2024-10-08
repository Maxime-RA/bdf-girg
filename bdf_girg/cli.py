import bdf_girg
import helper
from time import time
from bdfs.OneDimBDF import D
from bdfs.OuterMaxBDF import OuterMax
from bdfs.OuterMinBDF import OuterMin

help_string = 'Generate BDF-GIRG in expected linear time:\n' \
              '\t[-n anInt]    \t\t number of vertice                     \t default 1000\n' \
              '\t[-deg aFloat] \t\t average degree                        \t default 10\n' \
              '\t[-bdf aBDF]   \t\t boolean distance function             \t default min(0,max(1,2))\n' \
              '\t[-ple aFloat] \t\t power-law exponent                    \t default 2.5\n' \
              '\t[-comp]       \t\t ignore intersections                  \t default False\n' \
              '\t[-triv]       \t\t use trivial algorithm                 \t default False\n' \
              '\t[-simp]       \t\t use simplification algorithm          \t default False\n' \
              '\t[-pbdf]       \t\t stats about bdf                       \t default False\n' \
              '\t[-pmax]       \t\t stats about max-sets                  \t default False\n' \
              '\t[-wseed anInt]\t\t seed for weights                      \t default random\n' \
              '\t[-pseed anInt]\t\t seed for positions                    \t default random\n'\

using_parameter = 'boolean distance function : {}\n' \
                  'number of node            : {}\n' \
                  'average degree            : {}\n' \
                  'power-law exponent        : {}\n'


def print_stats():
    print("----- parameters -----------------------------")
    print(using_parameter.format(bdf, n, deg, ple))
    if pbdf:
        print("----- bdf stats ------------------------------")
        print("Computational depth  {} \t volumetric {}".format(bdf.get_depth_com(), bdf.get_depth_vol()))
        print("Computational length {} \t volumetric {}".format(bdf.get_length_com(), bdf.get_length_vol()))
        print("Min-max set: {}\n".format(bdf.get_min_max_form()))
    if pbdf and simp:
        s_bdf = bdf_girg.optimal_min_max_shortening(bdf)
        print(">>>>> simplified min-max set -------------------------")
        print("Depth  {} ".format(len(next(iter(s_bdf)))))
        print("Length {} ".format(len(s_bdf)))
        print("Min-max set: {}\n".format(s_bdf))


def print_girg_stats(edges, t):
    print('\t\t done in {:.3f}s'.format(time() - t))
    print("Generated   {} edges".format(len(edges)))
    av = (len(edges) * 2) / n
    print("Avg degree: {:.2f}  err:{:.2f}%\n".format(av, abs(((av / deg) - 1) * 100)))


def print_max_set_stats(pre_form, gen_form, post_form, t_pre, t_gen, t_pos, min_max_set):
    print(">>>> stats about bdf-girg generation -----")
    print("Total time for preprocessing:  {:.3f}s".format(t_gen - t_pre))
    print("Total time for computation:    {:.3f}s".format(t_pos - t_gen))
    print("Total time for postprocessing: {:.3f}s".format(time() - t_pos))
    print(">>>> stats about simplified max-set generation -----")
    for i in range(len(pre_form)):
        print("For set: {:<8}\n Time: pre. {:.3f}s  gen. {:.3f}s  post {:.3f}s \n Edges:  sampled {}. after check {}. added {}".format(
            str(min_max_set[i]),
            pre_form[i][2],
            gen_form[i][1],
            post_form[i][2],
            len(gen_form[i][0]),
            post_form[i][0],
            post_form[i][1]
        ))
    print("")


def generate_weights_positions():
    print("----- generate weights & position ------------")

    print('Generating positions...', end='', flush=True)
    t = time()
    gen_pos = bdf_girg.generate_positions(n, len(bdf.get_dimensions()), pseed)
    print('\t\t done in {:.3f}s'.format(time() - t))

    print('Generating weights...', end='', flush=True)
    t = time()
    gen_w = bdf_girg.generate_weights(n, ple, wseed)
    print('\t\t done in {:.3f}s\n'.format(time() - t))
    return gen_pos, gen_w


def estimate_threshold():
    print("----- estimating thresholds  ------------------")

    print('Estimate threshold constant...', end='', flush=True)
    t = time()
    constant = bdf_girg.estimate_threshold_constant(bdf, weights, deg, comp)
    print('\t done in {:.3f}s'.format(time() - t))

    if deg < bdf.get_depth_vol() * 4 * bdf.get_length_vol():
        t = time()
        print('Degree low, higher estimate...', end='', flush=True)
        const_high = bdf_girg.estimate_threshold_constant(bdf, weights, bdf.get_depth_vol() * 5 * bdf.get_length_vol(), comp)
        print('\t done in {:.3f}s'.format(time() - t))
    else:
        const_high = constant
    print("")
    return constant, const_high


def gen_trivial():
    print("----- generation using trivial algo -----------")

    print('Generating graph...', end='', flush=True)
    t = time()
    edges = bdf_girg.gen_bdf_edges_trivial(positions, weights, bdf, thr_con)
    print_girg_stats(edges, t)


def gen_simpl():
    print("-- generation using linear time algo ----------")
    print('Generating graph...', end='', flush=True)
    t_start = time()
    if not pmax:
        edges = bdf_girg.gen_bdf_girg(positions, weights, bdf, thr_con, thr_con_gen)
        print_girg_stats(set(edges), t_start)
        return
        
    MMS = list(bdf.get_min_max_form())
    t_pre = time()
    pre_form = bdf_girg.step_weight_adjust(positions, weights, bdf_girg.optimal_min_max_shortening(bdf), thr_con, bdf.get_depth_vol())
    
    t_gen = time()
    gen_form = bdf_girg.step_girg_gen(pre_form)
    
    t_pos = time()
    edges, pos_form = bdf_girg.step_girg_assemble(gen_form, positions, weights, MMS, thr_con, bdf.get_depth_vol())
    
    print_girg_stats(edges, t_start)
    print_max_set_stats(pre_form, gen_form, pos_form, t_pre, t_gen, t_pos, MMS)
    
    
   

if __name__ == '__main__':
    args = helper.parse_args()
    if not args:
        print(help_string)
        exit(0)

    # Set parameters
    n = int(args['n']) if 'n' in args else 1000
    deg = int(args['deg']) if 'deg' in args else 10
    bdf = helper.parse_bdf(args['bdf']) if 'bdf' in args else OuterMin(D(0), OuterMax(D(1), D(2)))
    ple = float(args['ple']) if 'ple' in args else 2.5
    comp = True if 'comp' in args else False
    triv = True if 'triv' in args else False
    simp = True if 'simp' in args else False
    pbdf = True if 'pbdf' in args else False
    pmax = True if 'pmax' in args else False
    wseed = int(args['wseed']) if 'wseed' in args else None
    pseed = int(args['pseed']) if 'pseed' in args else None

    assert set(bdf.get_dimensions()) == set(range(len(bdf.get_dimensions()))), "dimension should be from 0 to d"

    if not triv and not simp:
        print("Nothing to generate...")
        exit(0)

    print_stats()
    positions, weights = generate_weights_positions()
    thr_con, thr_con_gen = estimate_threshold()
    if triv:
        gen_trivial()
    if simp:
        gen_simpl()
