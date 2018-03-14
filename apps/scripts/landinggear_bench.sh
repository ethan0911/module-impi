#!/bin/bash

PREFIX=$1

$PREFIX/ospImplicitIsoSurfaceBench \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/cb.osp \
-iso 99000.0 \
-object ~feng/Desktop/ws/data/NASA_LandingGear_AMR/landingGear.obj \
-vp 16.103643 15.418482 0.265052 \
-vu -0.000000 -0.000000 -1.000000 \
-vi 16.387318 15.844191 0.332234 \
-sun 0.337 0.416 -0.605 \
-dis 0.554 1.000 -0.211 \
-scale 1 1 1 \
-translate 15.995 16 0.1 \
-valueRange 98280.0 99280.0 \
-frames 20 100 \
"${@:2}"
