#!/bin/bash

PREFIX=$1

$PREFIX/ospImplicitIsoSurfaceWidget \
/home/sci/feng/Desktop/ws/data/chombo/chombo_amr.osp \
-iso 0.7 \
-vp 24.684797 17.313093 -10.046009 \
-vu 0.000000 1.000000 0.000000 \
-vi 8.000000 7.999991 7.999997 \
-frames 50 200 \
"${@:2}"

