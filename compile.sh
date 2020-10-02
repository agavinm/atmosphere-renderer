#!/bin/bash

export CC=clang-10
export CXX=clang++-10
mkdir build
cd build
cmake -GNinja ..
ninja
cd ..
source setpath.sh
echo "$: mitsuba scenes/space/space.xml"
