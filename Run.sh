#!/usr/bin/sh

rm *.mp4
rm -rf data/*
mkdir data
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug
build/src/FourthOrderMicromagneticSolver