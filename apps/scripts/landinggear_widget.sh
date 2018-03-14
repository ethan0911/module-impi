#!/bin/bash

PREFIX=$1

$PREFIX/ospImplicitIsoSurfaceWidget \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/cb.osp \
-iso 99000.0 \
-object ~feng/Desktop/ws/data/NASA_LandingGear_AMR/landingGear.obj \
-vp 16.286070 16.446814 0.245150 \
-vu -0.000000 -0.000000 -1.000000 \
-vi 16.430407 16.157639 0.353916  \
-sun 0.337 0.416 -0.605 \
-dis 0.783 -1.00 -0.086 \
-scale 1 1 1 -translate 15.995 16 0.1 \
-frames 20 100 \
"${@:2}"
