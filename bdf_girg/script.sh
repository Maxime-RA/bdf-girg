#!/bin/bash
# Copy files to directory because of issues with dynamic linking
cd ..
make clean
make
cd girg_sampling/
cp _libgirgs_wrapper.so /home/maxime/Documents/GIRG/girg_python/bdf-girg/bdf_girg/
cp libgirgs.so.1 /home/maxime/Documents/GIRG/girg_python/bdf-girg/bdf_girg/
