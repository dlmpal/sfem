#!/bin/bash

mkdir -p build
rm -rf build/*

(cd ../../build; make)
(cd build; cmake -S .. -B .; make)

mkdir -p fields
rm fields/*

time mpiexec -np 4 build/laplaceSolver 3