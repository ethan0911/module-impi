#!/bin/bash
IMPI_AMR_DATA=landingGear \
IMPI_AMR_METHOD=octant \
./ospImplicitIsoSurfaceViewer \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/cb.osp \
/home/sci/feng/Desktop/ws/data/NASA_LandingGear_AMR/landingGear.obj \
-vp 16.286070 16.446814 0.245150 \
-vu -0.000000 -0.000000 -1.000000 \
-vi 16.430407 16.157639 0.353916 "$@"
