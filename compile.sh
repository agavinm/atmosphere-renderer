#!/bin/bash

export CC=clang-10
export CXX=clang++-10
mkdir build
cd build
cmake -DMTS_ENABLE_EMBREE=1 -GNinja ..
ninja
cd ..
source setpath.sh
