#!/bin/bash

set -e

cd ~

. ~/.bashrc

if [ -z "$(ls -A -- ~/build/)" ]; then
    echo Build directory is empty... Assuming we have to configure:
    echo
    cd ~/psi4/
    conda run --live-stream -n p4env conda/psi4-path-advisor.py cmake --objdir ../build --insist
    cat cache_p4env@build.cmake
    conda run --live-stream -n p4env cmake -S. -B ../build
    echo
    echo Configure complete!
    echo
else
    echo Build directory is NOT empty... Assuming the build is configured. Proceeding to build:
    echo
fi

cd ~/build/
conda run --live-stream -n p4env cmake --build . -j2

cd ~/work/
conda run --live-stream -n p4env snakemake -c1 --rerun-incomplete
