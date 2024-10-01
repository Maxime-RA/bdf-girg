# BDF-GIRG generation
An extension of the HyperGIRGs sampling algorithm for GIRGs enabling the sampling of the threshold model of BDF-GIRGs.
Arbitrary Boolean Distance Functions (BDF) are supported although higher dimensions (> 20) do not behave well with  threshold estimator due to numerical issues.  

A description of the algorithm can be found in [this](./thesis.pdf) bachlor thesis.
The paper describing the HyperGIRGs algorithm: [Efficiently Generating Geometric Inhomogeneous and Hyperbolic Random Graphs](https://arxiv.org/abs/1905.06706).

The C++ implementation of HyperGIRGs can be found [here](https://github.com/chistopher/girgs).
The Python wrapper we based our implementation on [here](https://github.com/gavento/girg-sampling).


## Install
First, make sure you have poetry and pybind11 installed. (`pip install poetry`, `pip install pybind11`).
Next, change the path in the Makefile (line 6) to where your pybind is installed.
Finally, in `/bdf-girg`, execute `./script` to compile everything.

## Usage
You can either use the command line interface with `python3 cli.py` or the module `bdf_girg` and `helper` to parse bdfs.
``` python
import bdf_girg, helper

n = 2 ** 15
degree = 10
bdf = helper.parse_bdf("min(0,max(1,2))"):
ple = 2.5

# Generate positions and weights
positions = bdf_girg.generate_positions(n, len(bdf.get_dimensions()))
weights = bdf_girg.generate_weights(n,ple)

# Estimate the threshold constant, not ignoring the intersections
thr_con = estimate_threshold_constant(bdf, weights, degree, False)

# Sample BDF-GIRG using the trivial O(n^2) algorithm 
edges_trivial = bdf_girg.gen_bdf_edges_trivial(positions, weights, bdf, thr_con)

# Sample edges in linear time
edges_linear = bdf_girg.gen_bdf_girg(positions, weights, bdf, thr_con)

# Get bdf shortening
MMS = list(bdf.get_min_max_form())
MMS_reduced = bdf_girg.optimal_min_max_shortening(bdf)

# Sample max-girgs seperatly to know for each by max-girg the number of 
# sampled and added edges.
pre_list = step_weight_adjust(positions, weights, MMS_reduced, thr_con, bdf.get_depth_vol())
gen_list = step_girg_gen(pre_list)
edges_stats, stats = bdf_girg.step_girg_assemble(gen_list, positions, weights, MMS, thr_con, bdf.get_depth_vol())


```
