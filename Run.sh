#!/usr/bin/sh

cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --config Debug
build/src/FourthOrderMicromagneticSolver