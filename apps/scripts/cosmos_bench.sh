#!/bin/bash

IMPI_AMR_METHOD=octant \
./ospImplicitIsoSurfaceBench \
/home/sci/feng/Desktop/ws/data/COSMOS_ospray/cosmos.osp \
-iso 0.0 \
-vp 814.398254 -54.506733 878.914001 \
-vu 1.000000 0.000000 0.000000 \
-vi 248.845215 267.241455 256.334686 \
-frames 50 200 \
"$@"
