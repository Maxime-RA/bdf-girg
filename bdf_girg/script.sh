#!/bin/bash
# Builds the wrapper and makes copies to correct directory due to issues with dynamic linking :/
cd ..
make clean
make
cd girg_sampling/
cp _libgirgs_wrapper.so ../bdf_girg/
cp libgirgs.so.1 ../bdf_girg/
